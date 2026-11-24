program cube_lfm
 use iso_fortran_env,only:real64
 use iso_c_binding,only:c_int,c_char,c_null_char
 use cube_event_global,only:global_write
 use cube_compact_output,only:cube_output
 use mpi_f08
 use lfm_cube_evolution
 use lfm_cube_yy_background,only:cube_yy_background_configure
 use lfm_cube_yy_mapping,only:yy_mapping_configure
 use lfm_cube_yy_boundary,only:yy_boundary_configure
 implicit none
 include 'yy_rotation_profile.inc'
 type(cube_output)::archive
 type(cube_model)::c
 character(1024)::map,path,arg,restart,mi_reference,mi_background,wind_input
 character(64)::fingerprint
 integer::st,nang,mode,maxsteps,s,unit,profile_unit,mi_unit,mi_written,dims(3),r
 integer::fine_meta(11)
 integer,allocatable::all_fine_meta(:,:)
 real(real64)::fine_times(16)
 real(real64),allocatable::all_fine_times(:,:)
 real(real64)::stop,dt,errors(3),stats(6),wall,next_output,cadence,started,cp_cadence,next_cp
 real(real64)::profile(11),setup_s,setup_mi_s,io_start,step_start,advance_before,diagnostic_start
 interface
  integer(c_int)function event_open(path)bind(C,name='cube_event_wind_open')
   import c_int,c_char
   character(kind=c_char),intent(in)::path(*)
  end function
 end interface
 call MPI_Init(st)
 ! Dedicated Earth64 port profile. Magnetic tilt is in the released harmonic
 ! coefficients; do not apply a second arbitrary dipole tilt rotation.
 call cube_yy_background_configure(6380000._real64,1.e-9_real64,1._real64, &
  yy_spin_omega_rad_s,yy_spin_phase0_rad,st)
 if(st/=0) error stop 'released YY background profile'
 call yy_mapping_configure(6380000._real64,6500000._real64,yy_spin_omega_rad_s,yy_spin_phase0_rad,st)
 if(st/=0) error stop 'released YY mapping profile'
 call yy_boundary_configure(6380000._real64,yy_spin_omega_rad_s,st)
 if(st/=0) error stop 'released YY wall profile' 
 call get_command_argument(1,map);call get_command_argument(2,path)
 call get_command_argument(3,arg);read(arg,*) nang
 call get_command_argument(4,arg);read(arg,*) mode
 call get_command_argument(5,arg);read(arg,*) stop
 call get_command_argument(6,arg);read(arg,*) maxsteps
 cadence=600
 if(command_argument_count()>=7) then
  call get_command_argument(7,arg);read(arg,*) cadence
 endif
 call get_command_argument(8,restart);call get_command_argument(9,fingerprint)
 cp_cadence=cadence
 if(command_argument_count()>=10) then
  call get_command_argument(10,arg);read(arg,*) cp_cadence
 endif
 call get_command_argument(11,mi_reference);call get_command_argument(12,mi_background)
 dims=[2,1,1]
 if(command_argument_count()>=13) then
  if(command_argument_count()<15) error stop 'MPI layout requires all three optional arguments13..15'
  do s=1,3
   call get_command_argument(12+s,arg);read(arg,*,iostat=st) dims(s)
   if(st/=0) error stop 'integer MPI layout arguments required'
  enddo
 endif
 if(cadence<=0.or.cp_cadence<=0) error stop 'positive output/checkpoint cadences required'
 if(len_trim(fingerprint)/=64) error stop "64-character frozen provenance fingerprint required"
 next_output=cadence
 started=MPI_Wtime()
 call get_command_argument(16,wind_input,status=st)
 if(st/=0.or.len_trim(wind_input)==0)error stop 'OMNI HDF5 path required at argument 16'
 st=event_open(trim(wind_input)//c_null_char)
 if(st/=0)error stop 'dated OMNI HDF5 load'
 call c%init(trim(map),nang,mode,trim(mi_reference),trim(mi_background),mpi_dims=dims)

 if(c%rank==0) write(*,'(a,i0,a,3(i0,1x),a,3(i0,1x),a,f12.6)') &
  'CUBE_RUNTIME_LAYOUT ranks=',c%ranks,' dims=',c%layout%dims,' local=',c%n,' Boris_CA_km_s=',cube_boris_c_km_s
 call archive%open(c,trim(path))
 if(len_trim(restart)>0.and.trim(restart)/="-") call archive%restore(c,trim(restart),fingerprint)
 if(c%rank==0) write(*,'(a,3(es24.16,1x))') &
  'CUBE_YY_ROTATION period_s omega_rad_s phase0_rad=', &
  yy_spin_period_s,yy_spin_omega_rad_s,yy_spin_phase0_rad
 setup_s=MPI_Wtime()-started
 setup_mi_s=c%phase_seconds(7);c%phase_seconds=0
 c%halo%profile_seconds=0;c%guard_boundary_seconds=0;c%face_guard_seconds=0;c%staggered_guard_seconds=0;c%outer_seconds=0
#ifdef LFM_CUBE_NATIVE_MI
 if(c%mi_enabled) then
  c%mi%gather_seconds=0;c%mi%sample_seconds=0;c%mi%solve_seconds=0;c%mi%wall_seconds=0
 endif
#endif
 ! Native/Jupiter policy: output and CP observe natural MHD crossings.
 ! Reconstruct the next diagnostic clocks with the same roundoff predicate
 ! on cold start and restart; do not introduce additional MHD endpoints.
 next_output=cadence
 do while(clock_reached(c%time,next_output))
  next_output=next_output+cadence
 enddo
 next_cp=cp_cadence
 do while(clock_reached(c%time,next_cp))
  next_cp=next_cp+cp_cadence
 enddo
 if(c%rank==0) then
  open(newunit=unit,file=trim(path)//'-steps.csv',status='new')
  write(unit,'(a)') 'step,time_s,dt_s,wall_s,rho_min,c_min,v_max,B_max,divB_relative,mass_g,'// &
   'CT_step,shared_face,cfl_dt_s,event_limited,maxdt_limited'
  open(newunit=profile_unit,file=trim(path)//'-profile.csv',status='new')
  write(profile_unit,'(a)') 'step,time_s,wall_s,setup_s,setup_MI_s,guards_s,CFL_s,hydro_s,EMF_s,update_s,'// &
   'audit_s,MI_s,advance_s,IO_s,diagnostics_s,cleanup_s'
#ifdef LFM_CUBE_NATIVE_MI
  if(c%mi_enabled) then
   open(newunit=mi_unit,file=trim(path)//'-mi-updates.csv',status='new')
   write(mi_unit,'(a)') 'step,time_s,updates,CPCP_N_kV,CPCP_S_kV,FAC_max_uA_m2,global_final_residual,'// &
    'gather_s,sample_s,solve_s,wall_voltage_s,thin_shell_doubles,'// &
    'global_residual_copy_N,global_residual_copy_S,refinements_copy_N,refinements_copy_S,budget_copy_N,budget_copy_S,converged_copy_N,converged_copy_S,trace_failures'
  endif
#endif
 endif
 mi_written=0;call write_mi_update()
 diagnostic_start=MPI_Wtime();call c%diagnose(stats);c%diagnostic_seconds=MPI_Wtime()-diagnostic_start
 io_start=MPI_Wtime();if(c%rank==0)call global_write(trim(path),c%step,c%time)
   call archive%frame(c)
 c%io_seconds=c%io_seconds+MPI_Wtime()-io_start
 if(c%rank==0) write(unit,'(i0,14(",",es24.16))') c%step,c%time,0._real64,0._real64,stats, &
  0._real64,0._real64,0._real64,0._real64,0._real64
 do s=1,maxsteps
  if(c%time>=stop) exit
  step_start=MPI_Wtime();advance_before=c%phase_seconds(8)
  call c%advance(stop,dt,errors)
  ! Literal native/Jupiter snap_event_clock rule, only for the run endpoint.
  if(abs(c%time-stop)<=64._real64*epsilon(c%time)*max(1._real64,abs(c%time))) c%time=stop
  c%cleanup_seconds=c%cleanup_seconds+MPI_Wtime()-step_start-(c%phase_seconds(8)-advance_before)
  call c%couple();call write_mi_update()
  diagnostic_start=MPI_Wtime();call c%diagnose(stats)
  c%diagnostic_seconds=c%diagnostic_seconds+MPI_Wtime()-diagnostic_start
  wall=MPI_Wtime()-started
  if(c%rank==0) then
   io_start=MPI_Wtime()
   write(unit,'(i0,14(",",es24.16))') c%step,c%time,dt,wall,stats,errors(1:2),c%cfl_candidate_dt, &
    merge(1._real64,0._real64,dt<min(6._real64,c%cfl_candidate_dt)), &
    merge(1._real64,0._real64,dt==6._real64.and.c%cfl_candidate_dt>=6._real64)
   flush(unit)
   if(mod(c%step,100)==0) write(*,'(a,i0,a,f12.4,a,es12.4,a,f10.2)') &
     'CUBE_STEP ',c%step,' t=',c%time,' dt=',dt,' wall=',wall
   c%io_seconds=c%io_seconds+MPI_Wtime()-io_start
  endif
  if(clock_reached(c%time,next_output).or.c%time>=stop.or.s==maxsteps) then
   io_start=MPI_Wtime()
   if(c%rank==0)call global_write(trim(path),c%step,c%time)
   call archive%frame(c)
   do while(clock_reached(c%time,next_output))
    next_output=next_output+cadence
   enddo
   c%io_seconds=c%io_seconds+MPI_Wtime()-io_start
  endif
  if(clock_reached(c%time,next_cp).or.c%time>=stop.or.s==maxsteps) then
   io_start=MPI_Wtime()
   call archive%checkpoint(c,fingerprint)
   do while(clock_reached(c%time,next_cp))
    next_cp=next_cp+cp_cadence
   enddo
   c%io_seconds=c%io_seconds+MPI_Wtime()-io_start
  endif
  if(mod(c%step,100)==0.or.c%time>=stop.or.s==maxsteps) then
   call c%profile(profile)
   if(c%rank==0) then
    write(profile_unit,'(i0,15(",",es24.16))') c%step,c%time,MPI_Wtime()-started,setup_s,setup_mi_s,profile
    flush(profile_unit)
   endif
  endif
 enddo
 if(c%rank==0) then
#ifdef LFM_CUBE_NATIVE_MI
  if(c%mi_enabled) close(mi_unit)
#endif
  close(unit);close(profile_unit);write(*,*) 'CUBE_RUN_END steps=',c%step,' time=',c%time,' target=',stop
 endif
 ! One final per-rank instrumentation dump; never part of checkpoint state.
 fine_times=[c%halo%profile_seconds,c%guard_boundary_seconds,c%face_guard_seconds, &
   c%staggered_guard_seconds,c%phase_seconds(1),c%phase_seconds(8),c%outer_seconds]
 fine_meta=[c%layout%face,c%layout%first,c%n,c%halo%targets,c%halo%donors,c%halo%received, &
   count(c%halo%terms(:c%halo%targets)==1)]
 allocate(all_fine_times(16,c%ranks),all_fine_meta(11,c%ranks))
 call MPI_Gather(fine_times,16,MPI_DOUBLE_PRECISION,all_fine_times,16,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,st)
 call MPI_Gather(fine_meta,11,MPI_INTEGER,all_fine_meta,11,MPI_INTEGER,0,MPI_COMM_WORLD,st)
 if(c%rank==0) then
  open(newunit=profile_unit,file=trim(path)//'-guards-profile.csv',status='new')
  write(profile_unit,'(a)') 'rank,face,first_i,first_j,first_k,local_i,local_j,local_k,'// &
   'targets,donors,received,one_donor_targets,'// &
   'primitive_pack_s,primitive_MPI_s,sample_unpack_s,interpolation_s,physical_BC_s,face_reconstruct_s,'// &
   'face_outer_BC_s,CT_face_install_s,staggered_pack_s,staggered_MPI_s,staggered_unpack_s,guards_total_s,advance_s,'// &
   'kaiju_outer_faces_s,kaiju_outer_EMF_s,kaiju_nudge_s'
  do r=1,c%ranks
   write(profile_unit,'(i0,11(",",i0),16(",",es24.16))') r-1,all_fine_meta(:,r),all_fine_times(:,r)
  enddo
  close(profile_unit)
 endif
 call archive%close(c)
 call MPI_Finalize(st)
contains
 subroutine write_mi_update()
#ifdef LFM_CUBE_NATIVE_MI
  real(real64)::tick
  if(.not.c%mi_enabled.or.c%rank/=0.or.c%mi_updates==mi_written) return
  tick=MPI_Wtime()
  write(mi_unit,'(i0,",",es24.16,",",i0,11(",",es24.16),7(",",i0))') c%step,c%time,c%mi_updates, &
   c%mi%cpcp,c%mi%facmax,c%mi%worst_residual,c%mi%gather_seconds,c%mi%sample_seconds, &
   c%mi%solve_seconds,c%mi%wall_seconds,real(6*4*c%layout%global(2)*c%layout%global(3)*8,real64), &
   c%mi%solve_SER,c%mi%solve_ITER,c%mi%solve_ITMAX,c%mi%solve_converged, &
   count(c%mi%mapped_trace_status/=0.and.c%mi%mapped_trace_status/=-999)
  flush(mi_unit);mi_written=c%mi_updates
  c%io_seconds=c%io_seconds+MPI_Wtime()-tick
#endif
 end subroutine
end program
