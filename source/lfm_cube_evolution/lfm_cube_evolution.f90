! Six-patch evolution adapter around native LFM flux, EMF, CFL and CT kernels.
module lfm_cube_evolution
 use iso_fortran_env,only:real64,real32,int32,int64
 use iso_c_binding,only:c_int,c_double,c_size_t
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 use mpi_f08
 use lfm_cube_topology
 use lfm_cube_mesh
 use lfm_cube_mpi
 use lfm_cube_native_geometry
 use lfm_cube_yy_background,only:cube_yy_background_context,cube_yy_background_ready
 use lfm_cube_yy_mapping,only:yy_corotation_edge
 use yy_hydro_reference_basis,only:yy_reference_coefficients,yy_reference_at
 use lfm_cube_yy_boundary,only:yy_boundary_ghost
 use lfm_cube_yy_magnetic_wall,only:yy_magnetic_wall_prepare,yy_magnetic_wall_apply
 use lfm_cube_yy_fac_source,only:yy_source_prepare,yy_source_refresh
#ifdef LFM_CUBE_NATIVE_MI
 use lfm_cube_ionosphere
 use lfm_rotating_mi_geometry,only:mi_edge_lfm
#endif
 implicit none
 real(real64),save,private::equivalent_cost_s(8)=0
 ! Derived scheduling only: all96 owners keep complete MHD work.
 type(MPI_Comm),save,private::rotation_work_comm=MPI_COMM_NULL
 logical,save,private::rotation_work_ready=.false.
 integer,save,private::rotation_work_rank=-1,rotation_work_ranks=0
 real(real64),save,private::rotation_work_time=-1
 real(real64),allocatable,save,private::rotation_edge_points(:,:),rotation_edge_values(:),rotation_local_values(:)
 integer,allocatable,save,private::rotation_counts(:),rotation_offsets(:)

 real(real64),allocatable,save::direct_cell_bank(:,:,:),direct_face_bank(:,:,:),direct_edge_bank(:,:,:),direct_quad_bank(:,:,:)
 real(real64),allocatable,save::direct_quad_values(:,:,:,:,:)
 integer(c_int),allocatable,save::direct_cell_active(:),direct_face_active(:),direct_edge_active(:),direct_quad_active(:)
 private
 public::cube_model,earth_wind,earth_tail_mix,earth_outer_primitive
 public::cube_boris_c_km_s,cube_staggered_patch_guards
 public::clock_reached,cube_magnetic_reference
 public::kaiju_wind_weight,kaiju_outer_state,boundary_cell_center
 public::boundary_cartesian_curl,kaiju_outer_emf,kaiju_nudge_outer
 ! Immutable geometry-derived cache; rebuilt on init/own restart, no state/type/ABI change.
 real(real64),allocatable,save::yy_background_pressure_correction(:,:,:,:,:)
 ! Private exact time-dependent HYDRO reference; stored physical CT/m0 state unchanged.
 ! MI computes the perturbation current about the same complete planetary CT reference.
 real(real64),allocatable,save::yy_hydro_full_quad(:,:,:,:,:),yy_hydro_cell_shift(:,:,:,:)
 real(real64),allocatable,save::yy_hydro_quad_coef(:,:,:,:,:,:),yy_hydro_cell_coef(:,:,:,:,:)
 ! Same native B0/B1 interface extended to the prescribed complete planetary field.
 real(real64),allocatable,save::yy_bfield_face_coef(:,:,:,:,:),yy_bfield_edge_coef(:,:,:,:,:,:)
 real(real64),allocatable,save::yy_planetary_initial_phi(:,:,:,:)
 real(real64),allocatable,save::yy_bfield_face_reference(:,:,:,:),yy_bfield_edge_reference(:,:,:,:,:)
 logical,allocatable,save::yy_zero_cell(:,:,:),yy_zero_face(:,:,:),yy_zero_edge(:,:,:,:),yy_zero_quad(:,:,:,:)
 include 'yy_rotation_profile.inc'
 integer(int64),save,private::plane_words_compared=0
 integer,parameter::ng=4
#ifndef LFM_CUBE_BORIS_C_KM_S
#define LFM_CUBE_BORIS_C_KM_S 1500
#endif
 real(real64),parameter::cube_boris_c_km_s=real(LFM_CUBE_BORIS_C_KM_S,real64)
 real(real64),parameter::gamma=5._real64/3,ca=cube_boris_c_km_s*1.e5_real64
 ! Immutable, contiguous per-plane arguments for the unchanged native kernels.
 ! One allocatable component is one complete plane: taking a 4-D slab would
 ! make component/channel slices strided and trigger hidden wrapper copies.
 type cube_hydro_plane_geometry
  real(real64),allocatable::volume(:,:),area(:,:),frame(:,:,:),background(:,:,:)
 end type
 type cube_bfield_plane_geometry
  real(real64),allocatable::xyz(:,:,:),area(:,:,:),normal(:,:,:,:),projection(:,:,:)
  real(real64),allocatable::length(:,:),frame(:,:,:),background(:,:,:)
 end type
 type cube_sweep_geometry
  type(cube_hydro_plane_geometry),allocatable::hydro(:)
  type(cube_bfield_plane_geometry),allocatable::bfield(:)
 end type
 type cube_model
  type(cube_meridian)::m
  type(cube_patch)::p
  type(cube_layout)::layout
  type(cube_halo_plan)::halo
  type(cube_shared_plan)::faces,edges,fluxes
  ! Reconstructed point ghosts remain untouched. Only HYDRO's BASE-time
  ! CLIP/PDM references at a shared patch face use its actual other cell.
  real(real64),allocatable::hydro_base_reference(:,:)
  integer::n(3),rank,ranks,step=0,mode=0
  real(real64)::time=0,dt_old=-1,rp=1,sw(8),moment=0
  real(real64)::reference_uniform(3)=0,cfl_candidate_dt=0
  real(real64)::wind_switch_s=-1,wind_bz_after=0,wind_reference_offset_s=0
  logical::diagnostic_trace=.false.
  ! Opt-in until geometry and full OFF/ON payload gates have passed.
  logical::sweep_cache_enabled=.false.,sweep_cache_ready=.false.
  integer(int64)::sweep_cache_bytes=0
  ! Allocatable container also avoids gfortran's incorrect finalizer for a
  ! fixed array of nested allocatable derived types (qualified lifecycle test).
  type(cube_sweep_geometry),allocatable::sweep(:)
  ! Exclusive elapsed categories; hardware clocks are never checkpoint state.
  ! guards/predictor, CFL, hydro, EMF, CT+primitive+BZZ, audit, MI, advance total.
  real(real64)::phase_seconds(8)=0,io_seconds=0,diagnostic_seconds=0,cleanup_seconds=0
  real(real64)::outer_seconds(3)=0
  real(real64)::guard_boundary_seconds=0,face_guard_seconds(3)=0,staggered_guard_seconds(3)=0
  logical::mi_enabled=.false.
  real(real64)::mi_cadence=10,mi_next=10
  integer::mi_updates=0
#ifdef LFM_CUBE_NATIVE_MI
  type(cube_ionosphere)::mi
#endif
  real(real64),allocatable::q(:,:,:,:),history(:,:,:,:),pred(:,:,:,:),phi(:,:,:,:),phi_history(:,:,:,:)
  real(real64),allocatable::background_flux(:,:,:,:),background_quad(:,:,:,:,:),background_cell(:,:,:,:)
  real(real64),allocatable::background_edge(:,:,:,:,:)
  real(real64),allocatable::xv(:,:,:,:),area(:,:,:,:,:),vol(:,:,:),correction(:,:,:,:,:)
  real(real64),allocatable::area_magnitude(:,:,:,:) ! immutable exact norm2 of actual face vectors
  ! Static geometry only; not evolved or checkpoint state. Kaiju boundary adapter.
  real(real64),allocatable::outer_area(:,:,:,:,:),outer_xyz(:,:,:,:,:)
  real(real64),allocatable::outer_inverse(:,:,:,:,:),outer_curl_inverse(:,:,:,:,:)
  ! Immutable boundary geometry. Evaluate with the original expressions once;
  ! do not cache evolved fields, the convected wind front, dt or diffusion.
  real(real64),allocatable::outer_weight(:,:,:,:),outer_magnitude(:,:,:,:)
  real(real64),allocatable::outer_cell_xyz(:,:,:,:),outer_cell_weight(:,:,:)
  real(real64),allocatable::outer_normal(:,:,:)
  real(real64),allocatable::outer_edge_xyz(:,:,:,:,:),outer_edge_vector(:,:,:,:,:)
  real(real64),allocatable::outer_edge_weight(:,:,:,:),outer_edge_dl(:,:,:,:)
  real(real64)::outer_wind_speed=0
 contains
  procedure::configure_wind=>model_configure_wind
  procedure::init=>model_init
  procedure::prepare_sweeps=>model_prepare_sweeps
  procedure::advance=>model_advance
  procedure::guards=>model_guards
  procedure::face_guards=>model_face_guards
  procedure::diagnose=>model_diagnose
  procedure::dump=>model_dump
  procedure::checkpoint=>model_checkpoint
  procedure::restore=>model_restore
  procedure::profile=>model_profile
  procedure::couple=>model_couple
  procedure::dump_mi=>model_dump_mi
 end type
 interface
  integer(c_int) function pack_hydro_plane(n,g,d,k,q,p,s,qb,qp) bind(C,name='cube_yy_pack_hydro_plane')
   import c_int,c_double
   integer(c_int),intent(in)::n(3)
   integer(c_int),value::g,d,k
   real(c_double),intent(in)::q(*),p(*),s(*)
   real(c_double),intent(out)::qb(*),qp(*)
  end function
  integer(c_int) function pack_bfield_plane(n,g,d,k,q,p,bf,bp,ref,area,mag,u,emf,qb,qp,bb,pp) bind(C,name='cube_yy_pack_bfield_plane')
   import c_int,c_double
   integer(c_int),intent(in)::n(3)
   integer(c_int),value::g,d,k
   real(c_double),intent(in)::q(*),p(*),bf(*),bp(*),ref(*),area(*),mag(*),u(3),emf(3)
   real(c_double),intent(out)::qb(*),qp(*),bb(*),pp(*)
  end function
  integer(c_int) function distributed_wall_prepare(comm) bind(C,name='cube_yy_distributed_wall_prepare')
   import c_int
   integer(c_int),value::comm
  end function
  integer(c_int) function distributed_wall_prefetch(comm,time,dt_old) bind(C,name='cube_yy_distributed_wall_prefetch')
   import c_int,c_double
   integer(c_int),value::comm
   real(c_double),value::time,dt_old
  end function
  integer(c_int) function distributed_wall_check(comm,time) bind(C,name='cube_yy_distributed_wall_check')
   import c_int,c_double
   integer(c_int),value::comm
   real(c_double),value::time
  end function
  subroutine distributed_wall_release() bind(C,name='cube_yy_distributed_wall_release')
  end subroutine
  integer(c_int) function direct_reference_eval(count,bank,active,w,r) bind(C,name='cube_yy_direct_reference_evaluate')
   import c_int,c_double
   integer(c_int),value::count
   real(c_double),intent(in)::bank(*),w(*)
   integer(c_int),intent(in)::active(*)
   real(c_double),intent(out)::r(*)
  end function
  subroutine cube_yy_background_pressure_delta(v,area,delta)
   import real64
   real(real64),intent(in)::v(3,4),area(3)
   real(real64),intent(out)::delta(3)
  end subroutine
  integer(c_int) function yy_flux_trap(ni,nj,nk,f) bind(C,name='cube_yy_inner_fluid_flux_boundary')
   import c_int,c_double
   integer(c_int),value::ni,nj,nk
   real(c_double),intent(inout)::f(*)
  end function
  subroutine yy_flux_trap_statistics(calls,faces,massmax,energymax) bind(C,name='cube_yy_inner_flux_statistics')
   import c_size_t,c_double
   integer(c_size_t),intent(out)::calls,faces
   real(c_double),intent(out)::massmax,energymax
  end subroutine
  subroutine lfm_cube_hydro_flux(ni,nj,g,dt,gg,cc,pdm,base,pred,vol,area,frame,bn,bg,flux,acc,status)
   import real64
   integer,intent(in)::ni,nj,g
   real(real64),intent(in)::dt,gg,cc,pdm
   real(real64),intent(in)::base(1-g:ni+g,1-g:nj+g,8),pred(1-g:ni+g,1-g:nj+g,8)
   real(real64),intent(in)::vol(1-g:ni+g,1-g:nj+g),area(ni+1,nj),frame(ni+1,nj,9),bn(ni+1,nj),bg(ni+1,nj,8)
   real(real64),intent(out)::flux(ni+1,nj,8),acc(ni,nj,8)
   integer,intent(out)::status
  end subroutine
  subroutine lfm_cube_bfield(ni,nj,g,base,pred,xyz,edge,frame,area,normal,bface,bpred,bg,emf,status)
   import real64
   integer,intent(in)::ni,nj,g
   real(real64),intent(in)::base(1-g:ni+g,1-g:nj+g,8),pred(1-g:ni+g,1-g:nj+g,8)
   real(real64),intent(in)::xyz(1-g:ni+g,1-g:nj+g,3),edge(ni+1,nj+1),frame(ni+1,nj+1,9)
   real(real64),intent(in)::area(1-g:ni+g,1-g:nj+g,2),normal(1-g:ni+g,1-g:nj+g,3,2)
   real(real64),intent(in)::bface(1-g:ni+g,1-g:nj+g,2),bpred(1-g:ni+g,1-g:nj+g,2),bg(ni+1,nj+1,3)
   real(real64),intent(out)::emf(ni+1,nj+1)
   integer,intent(out)::status
  end subroutine
 end interface
contains
 pure logical function clock_reached(t,event) result(reached)
  real(real64),intent(in)::t,event
  ! Same scheduler convention as the unchanged original MPI Earth driver.
  reached=t>=event-64._real64*epsilon(t)*max(1._real64,abs(t))
 end function
 function earth_wind(sw,time,x,switch_s,bz_after,reference_offset_s) result(q)
  real(real64),intent(in)::sw(8),time,x(3)
  real(real64),optional,intent(in)::switch_s,bz_after,reference_offset_s
  real(real64)::q(8),tm,schedule_tm
  ! Original IBOUND_FOR's time4, constant solar-wind table and convected
  ! IMF front. Its single-precision time conversion is retained explicitly.
  q=sw;tm=real(real(time,real32),real64)-dot_product(sw(2:4),x)/sum(sw(2:4)**2)
  if(nint(tm)<=0) q(7:8)=0
  if(present(switch_s).neqv.present(bz_after)) error stop 'paired wind schedule arguments required'
  if(present(reference_offset_s).and..not.present(switch_s)) error stop "wind offset requires schedule"
  if(present(switch_s).and.present(bz_after)) then
   schedule_tm=tm
   if(present(reference_offset_s)) schedule_tm=tm+reference_offset_s
   ! Original IBOUND_FOR's one-second NINT table lookup, with the native
   ! schedule table rule BZ(I)=after when I>=switch. The optional shift
   ! is retained only for explicit legacy upstream-reference inputs. The
   ! qualified cold IMF front remains in its original X=0 convention.
   if(switch_s>0.and.nint(tm)>0) then
    if(real(nint(schedule_tm),real64)>=switch_s) q(8)=bz_after
   endif
  endif
 end function
 real(real64) function earth_tail_mix(x) result(w)
  real(real64),intent(in)::x(3)
  real(real64)::cone
  w=0
  if(x(1)>=0) return
  cone=sqrt((sum(x(2:)**2)+1.e-20_real64)/(sum(x*x)+1.e-20_real64))
  w=.5_real64*(1-tanh((cone-.3_real64)*12._real64))
 end function
 function earth_outer_primitive(sw,time,x,donor,layer,switch_s,bz_after,reference_offset_s) result(q)
  real(real64),intent(in)::sw(8),time,x(3),donor(8)
  real(real64),optional,intent(in)::switch_s,bz_after,reference_offset_s
  integer,intent(in)::layer
  real(real64)::q(8),cone,drop
  integer::i
  cone=(sum(x(2:)**2)+1.e-20_real64)/(sum(x*x)+1.e-20_real64)
  if(x(1)<0.and.cone<.09_real64) then
   drop=.9_real64
   do i=2,layer
    drop=.9_real64*drop
   enddo
   q=donor;q(1)=q(1)*drop
  else
   q=earth_wind(sw,time,x,switch_s,bz_after,reference_offset_s)
  endif
 end function
 subroutine demand(ok,msg)
  logical,intent(in)::ok
  character(*),intent(in)::msg
  integer::ierr
  if(ok) return
  write(*,*) 'CUBE_EVOLUTION_FAIL ',msg
  call MPI_Abort(MPI_COMM_WORLD,1,ierr)
 end subroutine
 subroutine frame_from_normal(normal,t)
  real(real64),intent(in)::normal(3)
  real(real64),intent(out)::t(9)
  real(real64)::seed(3)
  integer::a
  t(1:3)=normal/norm2(normal);a=minloc(abs(t(1:3)),dim=1);seed=0;seed(a)=1
  t(4:6)=seed-dot_product(seed,t(1:3))*t(1:3);t(4:6)=t(4:6)/norm2(t(4:6))
  t(7:9)=cube_cross(t(1:3),t(4:6))
 end subroutine
 subroutine extended_point(m,face,n,g,x)
  type(cube_meridian),intent(in)::m
  integer,intent(in)::face,n,g(3)
  real(real64),intent(out)::x(3)
  real(real64)::near(3),r
  integer::v(3),s
  if(m%shell%ns>0) then
   call m%point(face,n,g,x,s);call demand(s==0,'logical shell extended geometry');return
  endif
  v=g;v(1)=max(0,min(m%ni,g(1)))
  call m%point(face,n,v,x,s);call demand(s==0,'extended geometry ray')
  if(g(1)==v(1)) return
  v(1)=merge(1,m%ni-1,g(1)<0)
  call m%point(face,n,v,near,s);call demand(s==0,'extended geometry adjacent ray')
  r=norm2(near)/norm2(x)
  if(g(1)<0) then
   x=x*r**g(1)
  else
   x=x*(1/r)**(g(1)-m%ni)
  endif
 end subroutine
 ! Optional immutable forcing input; derived phase has no restart history.
 ! The checkpoint's run fingerprint includes this exact input file hash.
 subroutine model_configure_wind(c,path)
  class(cube_model),intent(inout)::c
  character(*),intent(in)::path
  character(80)::magic,extra
  real(real64)::value(2),minimum(2),maximum(2),local_offset,speed2
  integer::u,st,i,j,k
  logical::earth_reference
  call demand(c%mode>=2,'wind schedule is an Earth boundary configuration')
  open(newunit=u,file=path,status='old',action='read',iostat=st)
  call demand(st==0,'open frozen wind schedule')
  read(u,'(a)',iostat=st) magic
  call demand(st==0,'read wind schedule version')
  earth_reference=trim(magic)=='LFM_CUBE_IMF_STEP_EARTH_V1'
  call demand(earth_reference.or.trim(magic)=='LFM_CUBE_IMF_STEP_UPSTREAM_V1','wind schedule version')
  read(u,*,iostat=st) value
  call demand(st==0.and.all(ieee_is_finite(value)),'finite wind schedule values')
  call demand(value(1)>0.and.value(1)<real(huge(0)-2,real64),'positive representable switch time')
  call demand(value(1)==real(nint(value(1)),real64),'native wind table uses integer seconds')
  read(u,'(a)',iostat=st) extra;close(u)
  call demand(st<0,'unexpected wind schedule suffix')
  call MPI_Allreduce(value,minimum,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,st)
  call demand(st==MPI_SUCCESS,'wind schedule minimum collective')
  call MPI_Allreduce(value,maximum,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,st)
  call demand(st==MPI_SUCCESS.and.all(minimum==maximum),'identical wind schedule on all ranks')
  ! Native IBOUND_FOR: table epochs refer to x=0. Its existing convected
  ! tm=t-v.x/v^2 already injects the change early on the positive-X boundary.
  call demand(allocated(c%p%x),'active immutable geometry required for convected clock')
  speed2=sum(c%sw(2:4)**2);call demand(speed2>0,'nonzero convected wind speed')
  c%wind_reference_offset_s=0
  if(.not.earth_reference) then
   local_offset=huge(local_offset)
   do k=lbound(c%p%x,4),ubound(c%p%x,4)
    do j=lbound(c%p%x,3),ubound(c%p%x,3)
     do i=lbound(c%p%x,2),ubound(c%p%x,2)
      local_offset=min(local_offset,dot_product(c%sw(2:4),c%p%x(:,i,j,k))/speed2)
     enddo
    enddo
   enddo
   call MPI_Allreduce(local_offset,c%wind_reference_offset_s,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,st)
   call demand(st==MPI_SUCCESS.and.ieee_is_finite(c%wind_reference_offset_s),'common upstream clock from actual mesh')
  endif
  c%wind_switch_s=value(1);c%wind_bz_after=value(2)
  if(c%rank==0) then
   write(*,'(a,a)') 'CUBE_WIND_TIME_REFERENCE ',trim(magic)
   write(*,'(a,4es25.16)') &
    'CUBE_WIND_SCHEDULE table_s,BzBefore_G,BzAfter_G,convected_ref_offset_s: ', &
    c%wind_switch_s,c%sw(8),c%wind_bz_after,c%wind_reference_offset_s
  endif
 end subroutine
 subroutine model_init(c,path,nang,mode,mi_reference,mi_background,mpi_dims)
  class(cube_model),intent(out)::c
  character(*),intent(in)::path
  integer,intent(in)::nang,mode
  character(*),optional,intent(in)::mi_reference,mi_background
  integer,optional,intent(in)::mpi_dims(3)
  integer::st,i,j,k,d,a,b,h,n(3),idx(3),vtx(3,4),z,decomposition(3)
  character(32)::trace
  real(real64),allocatable::ref_b0(:,:,:,:)
  real(real64)::v(3,0:1,0:1,0:1),delta(3,3),v4(3,4),mag(3),x(3),aa(3),bb(3),t,condition
  call MPI_Comm_rank(MPI_COMM_WORLD,c%rank,st);call MPI_Comm_size(MPI_COMM_WORLD,c%ranks,st)
  call get_environment_variable('CUBE_DIAGNOSTICS',trace,status=st)
  c%diagnostic_trace=st==0.and.trim(trace)=='1'
  call get_environment_variable('CUBE_SWEEP_CACHE',trace,status=st)
  c%sweep_cache_enabled=st==0.and.trim(trace)=='1'
  if(c%rank==0) write(*,*) 'INIT_START'
  call c%m%load(path,st);call demand(st==0,'input meridian')
  if(c%rank==0) write(*,*) 'INIT_MAP_LOAD_PASS'
  call demand(c%m%ni<=128.and.c%ranks<=128,'runtime limits radial<=128 and MPI ranks<=128')
  call demand(cube_boris_c_km_s>0.and.ieee_is_finite(cube_boris_c_km_s),'positive compiled Boris speed')
  decomposition=[2,1,1]
  if(present(mpi_dims)) decomposition=mpi_dims
  c%mode=mode
  call demand(mode>=0.and.mode<=3,'supported model mode0/1/2/3')
  c%mi_enabled=mode==3
#ifndef LFM_CUBE_NATIVE_MI
  call demand(.not.c%mi_enabled,'mode3 requires linked original native ionosphere')
#endif
  if(mode==1) then
   c%rp=7.1492e9_real64;c%moment=-4.3_real64*c%rp**3
   c%sw=[3.34524384e-25_real64,-4.e7_real64,0._real64,0._real64,4.e6_real64, &
      0._real64,0._real64,5.e-5_real64]
   call c%m%scale(c%rp)
  else if(mode>=2) then
   c%rp=6.38e8_real64;c%moment=.30_real64*c%rp**3
   c%sw=[8.3631096e-24_real64,-4.e7_real64,0._real64,0._real64,4.e6_real64, &
      0._real64,0._real64,5.e-5_real64]
   call c%m%scale(c%rp)
  else
   c%rp=1;c%moment=0
   c%sw=[1._real64,100._real64,30._real64,20._real64,10._real64,.2_real64,.4_real64,.1_real64]
  endif
  c%reference_uniform=c%sw(6:8)
  if(mode>=2) c%reference_uniform=0
  call c%layout%init([c%m%ni,nang,nang],decomposition,c%rank,c%ranks,st)
  call demand(st==0,'six-patch regular MPI layout')
  c%n=c%layout%local;n=c%n
  call demand(all(n>=ng+1),'every local extent must be at least ghost width plus1')
  if(c%rank==0) write(*,*) 'INIT_LAYOUT_PASS ',n
  call c%p%init_mapped(c%layout,c%m,st);call demand(st==0,'active mapped geometry')
  if(c%rank==0) write(*,*) 'INIT_ACTIVE_MESH_PASS'
  call native_load(c%p)
  if(c%rank==0) write(*,*) 'INIT_NATIVE_METRIC_PASS'
  call cube_bzz_calibration(c%p,c%correction,condition,st)
  call demand(st==0,'native BZZ calibration')
  allocate(c%xv(3,-ng:n(1)+ng,-ng:n(2)+ng,-ng:n(3)+ng))
  do k=-ng,n(3)+ng
   do j=-ng,n(2)+ng
    do i=-ng,n(1)+ng
     call extended_point(c%m,c%layout%face,nang,c%layout%first+[i,j,k],c%xv(:,i,j,k))
    enddo
   enddo
  enddo
  allocate(c%q(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng,8));allocate(c%history,mold=c%q)
  allocate(c%pred,mold=c%q)
  allocate(c%area(3,1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng,3))
  allocate(c%area_magnitude(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng,3))
  allocate(c%vol(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng))
  do k=1-ng,n(3)+ng
   do j=1-ng,n(2)+ng
    do i=1-ng,n(1)+ng
     do h=0,7
      v(:,mod(h,2),mod(h/2,2),h/4)=c%xv(:,i-1+mod(h,2),j-1+mod(h/2,2),k-1+h/4)
     enddo
     delta(:,1)=sum(sum(v(:,1,:,:)-v(:,0,:,:),dim=3),dim=2)/4
     delta(:,2)=sum(sum(v(:,:,1,:)-v(:,:,0,:),dim=3),dim=2)/4
     delta(:,3)=sum(sum(v(:,:,:,1)-v(:,:,:,0),dim=3),dim=2)/4
     c%vol(i,j,k)=abs(dot_product(delta(:,1),cube_cross(delta(:,2),delta(:,3))))
     c%q(i,j,k,:)=c%sw
     if(mode>=1) then
      x=sum(sum(sum(v,dim=4),dim=3),dim=2)/8
      if(x(1)<200*c%rp.or.mode>=2) then
       c%q(i,j,k,1)=2.25e-25_real64;c%q(i,j,k,2:4)=0;c%q(i,j,k,5)=1.5e6_real64
      endif
     endif
     do d=1,3
      idx=[i,j,k];call get_face(c,d,idx,v4);c%area(:,i,j,k,d)=cube_quad(v4)
      c%area_magnitude(i,j,k,d)=norm2(c%area(:,i,j,k,d))
     enddo
    enddo
   enddo
  enddo
  call demand(minval(c%vol)>0,'positive active and guard volumes')
  call demand(all(ieee_is_finite(c%area_magnitude)).and.minval(c%area_magnitude)>0,'positive fixed face norms')
  if(c%mode>=1.and.c%layout%coord(1)==0) then
   call yy_magnetic_wall_prepare(ng,n,c%xv,c%area,c%area_magnitude,st)
   call demand(st==0,'released rotating magnetic wall fixed original Cube geometry')
  endif
  if(c%mode>=1) then
   call yy_source_prepare(c%m,c%layout,ng,c%xv,.01_real64*c%rp,MPI_COMM_WORLD,st)
   call demand(st==0,'released current immutable original Cube geometry')
  endif
  allocate(c%phi(0:n(1),0:n(2),0:n(3),3));allocate(c%phi_history,mold=c%phi)
  c%phi=0;mag=c%reference_uniform
  do d=1,3
   do k=0,n(3)
    do j=0,n(2)
     do i=0,n(1)
      idx=[i,j,k]
      if(any(pack(idx,[d/=1,d/=2,d/=3])==0)) cycle
      idx(d)=idx(d)+1;c%phi(i,j,k,d)=dot_product(mag,c%area(:,idx(1),idx(2),idx(3),d))
     enddo
    enddo
   enddo
  enddo
  call c%faces%init(c%p,4,MPI_COMM_WORLD,st);call demand(st==0,'face ownership')
  call c%edges%init(c%p,2,MPI_COMM_WORLD,st);call demand(st==0,'edge ownership')
  call c%fluxes%init(c%p,4,MPI_COMM_WORLD,st,components=8,reconcile=.true.)
  call demand(st==0,'shared hydro flux topology plan')
  call c%fluxes%prepare_cell_pairs(c%n,st)
  call demand(st==0,'physical adjacent base-reference plan')
  allocate(c%hydro_base_reference(8,c%fluxes%entries))
  call demand(st==0,'eight-channel fluid seam reconciliation')
  if(c%rank==0) write(*,*) 'INIT_SHARED_PLANS_PASS'
  call c%halo%init(c%layout,c%m,ng,8,MPI_COMM_WORLD,st);call demand(st==0,'physical cell guards')
  allocate(c%background_flux,mold=c%phi);c%background_flux=0
  allocate(c%background_quad(8,0:n(1),0:n(2),0:n(3),3));c%background_quad=0
  call demand(.not.allocated(yy_background_pressure_correction),'one immutable background pressure cache per model')
  allocate(yy_background_pressure_correction(3,0:n(1),0:n(2),0:n(3),3));yy_background_pressure_correction=0
  allocate(c%background_cell(n(1),n(2),n(3),3));c%background_cell=0
  allocate(c%background_edge(3,0:n(1),0:n(2),0:n(3),3));c%background_edge=0
  call lfm_cube_background_setup(c%moment)
  if(mode>=2) call lfm_cube_earth_background(c%moment,c%rp)
  if(mode>=1) then
   call add_initial_dipole(c)
   call background_quadratures(c)
  endif
  call prepare_hydro_time_reference(c)
  call export_static_background(c)
  call c%faces%sync(c%phi,st);call demand(st==0,'initial shared face flux')
  call native_recover(c%p,c%phi-c%background_flux,c%correction,c%q(1:n(1),1:n(2),1:n(3),6:8))
  if(c%rank==0) write(*,*) 'INIT_RECOVERY_PASS'
  call prepare_outer_boundary(c)
  call c%guards(c%q);c%history=c%q;c%phi_history=c%phi
  if(c%sweep_cache_enabled) call c%prepare_sweeps()
#ifdef LFM_CUBE_NATIVE_MI
  if(c%mi_enabled) then
   call demand(present(mi_reference),'native MI reference geometry required')
   ! Released YY coupling takes curl of the advected residual. Do not add
   ! the old dipole reference field to the sampled IGRF residual packet.
   call c%mi%init(c%m,c%layout,trim(mi_reference),MPI_COMM_WORLD,st)
   call demand(st==0,'native ionosphere fixed reference/cached cube sampler initialization')
   c%mi_next=0;call c%couple()
  endif
#endif
  call prepare_rotation_work(c)
 end subroutine
 subroutine prepare_rotation_work(c)
  class(cube_model),intent(in)::c
  integer::colour,key,st,i,j,k,l,nedge,local,global,owned
  real(real64),parameter::epochs(8)=[0._real64,.173456789_real64,7200.0050561318985_real64, &
   10800._real64,21600.104527889976_real64,32400._real64,43199.823456789_real64,43200._real64]
  real(real64)::original
  integer::epoch
  if(.not.c%mi_enabled.or.c%layout%dims(1)/=4)return
  colour=c%layout%face+6*(c%layout%coord(2)+c%layout%dims(2)*c%layout%coord(3));key=c%layout%coord(1)
  call MPI_Comm_split(MPI_COMM_WORLD,colour,key,rotation_work_comm,st)
  call demand(st==MPI_SUCCESS,'private four-owner rotation communicator')
  call MPI_Comm_rank(rotation_work_comm,rotation_work_rank,st)
  call demand(st==MPI_SUCCESS.and.rotation_work_rank==key,'radial rotation worker owner identity')
  call MPI_Comm_size(rotation_work_comm,rotation_work_ranks,st)
  call demand(st==MPI_SUCCESS.and.rotation_work_ranks==4,'four existing radial MHD owners')
  nedge=c%n(2)*(c%n(3)+1)+(c%n(2)+1)*c%n(3)
  allocate(rotation_edge_points(6,nedge),rotation_edge_values(nedge),rotation_counts(4),rotation_offsets(4),stat=local)
  call MPI_Allreduce(local,global,1,MPI_INTEGER,MPI_MAX,rotation_work_comm,st)
  call demand(st==MPI_SUCCESS.and.global==0,'rotation derived-buffer allocation')
  do i=1,4
   rotation_offsets(i)=nedge*(i-1)/4;rotation_counts(i)=nedge*i/4-rotation_offsets(i)
  enddo
  owned=rotation_counts(rotation_work_rank+1);allocate(rotation_local_values(owned),stat=local)
  call MPI_Allreduce(local,global,1,MPI_INTEGER,MPI_MAX,rotation_work_comm,st)
  call demand(st==MPI_SUCCESS.and.global==0,'rotation local-buffer allocation')
  if(rotation_work_rank==0)then
   l=0
   do k=0,c%n(3);do j=1,c%n(2)
    l=l+1;rotation_edge_points(1:3,l)=.01_real64*c%xv(:,0,j-1,k)
    rotation_edge_points(4:6,l)=.01_real64*c%xv(:,0,j,k)
   enddo;enddo
   do k=1,c%n(3);do j=0,c%n(2)
    l=l+1;rotation_edge_points(1:3,l)=.01_real64*c%xv(:,0,j,k-1)
    rotation_edge_points(4:6,l)=.01_real64*c%xv(:,0,j,k)
   enddo;enddo
   call demand(l==nedge,'same original rotation edge order')
  endif
  call MPI_Bcast(rotation_edge_points,size(rotation_edge_points),MPI_DOUBLE_PRECISION,0,rotation_work_comm,st)
  call demand(st==MPI_SUCCESS.and.all(ieee_is_finite(rotation_edge_points)),'exact original inner-edge coordinate broadcast')
  st=distributed_wall_prepare(int(rotation_work_comm%MPI_VAL,c_int))
  call demand(st==0,'original magnetic-wall geometry and derived cache worker broadcast')
  rotation_work_ready=.true.
  ! New execution-owner qualification only. No physical clock/default changes.
  do epoch=1,size(epochs)
   st=distributed_wall_prefetch(int(rotation_work_comm%MPI_VAL,c_int),epochs(epoch),-1._real64)
   call demand(st==0,'original magnetic-wall worker epoch status')
   st=distributed_wall_check(int(rotation_work_comm%MPI_VAL,c_int),epochs(epoch))
   call demand(st==0,'original magnetic-wall all-point all-worker WORD equality')
   call refresh_rotation_work(epochs(epoch))
   local=0
   if(rotation_work_rank==0)then
    do l=1,nedge
     call yy_corotation_edge(rotation_edge_points(1:3,l),rotation_edge_points(4:6,l),epochs(epoch),original,st)
     if(st/=0)then
      local=1;exit
     endif
     if(transfer(original,0_int64)/=transfer(rotation_edge_values(l),0_int64))then
      write(*,*)'ROTATION_WORKER_WORD_MISMATCH ',c%rank,epoch,l,original,rotation_edge_values(l)
      local=2;exit
     endif
    enddo
   endif
   call MPI_Bcast(local,1,MPI_INTEGER,0,rotation_work_comm,st)
   call demand(st==MPI_SUCCESS.and.local==0,'original corotation function all-edge all-worker WORD equality')
  enddo
  rotation_work_time=-1
  if(rotation_work_rank==0)write(*,*)'ORIGINAL_COROTATION_DISTRIBUTED_WORD_PASS ',c%rank,nedge*size(epochs)
 end subroutine
 subroutine refresh_rotation_work(time)
  real(real64),intent(in)::time
  integer::first,last,l,j,local,global,st
  if(.not.rotation_work_ready)return
  first=rotation_offsets(rotation_work_rank+1)+1;last=first+rotation_counts(rotation_work_rank+1)-1
  local=0;j=0
  do l=first,last
   j=j+1
   call yy_corotation_edge(rotation_edge_points(1:3,l),rotation_edge_points(4:6,l),time,rotation_local_values(j),st)
   if(st/=0)then
    local=1;exit
   endif
  enddo
  call MPI_Allreduce(local,global,1,MPI_INTEGER,MPI_MAX,rotation_work_comm,st)
  call demand(st==MPI_SUCCESS.and.global==0,'original rotating-edge function worker status')
  call MPI_Gatherv(rotation_local_values,size(rotation_local_values),MPI_DOUBLE_PRECISION,rotation_edge_values, &
   rotation_counts,rotation_offsets,MPI_DOUBLE_PRECISION,0,rotation_work_comm,st)
  call demand(st==MPI_SUCCESS,'exact rotating-edge result gather')
  rotation_work_time=time
 end subroutine
 subroutine release_rotation_work()
  integer::st
  if(.not.rotation_work_ready)return
  call distributed_wall_release()
  call MPI_Comm_free(rotation_work_comm,st)
  call demand(st==MPI_SUCCESS,'private rotating-edge communicator release')
  rotation_work_ready=.false.;rotation_work_time=-1
  deallocate(rotation_edge_points,rotation_edge_values,rotation_local_values,rotation_counts,rotation_offsets)
 end subroutine
 subroutine model_prepare_sweeps(c)
  class(cube_model),intent(inout)::c
  real(real64)::normal(3),frame(9),edge(3),a0(3),b0(3),pm(3)
  integer::n(3),d,a,baxis,kk,ii,jj,cvar,i,idx(3),ed(3),low(3),high(3)
  if(c%sweep_cache_ready) return
  n=c%n;c%sweep_cache_bytes=0
  if(.not.allocated(c%sweep)) allocate(c%sweep(3))
  call demand(size(c%sweep)==3,'three immutable sweep directions')
  do d=1,3
   a=mod(d,3)+1;baxis=mod(d+1,3)+1
   allocate(c%sweep(d)%hydro(n(baxis)),c%sweep(d)%bfield(n(d)))
   do kk=1,n(baxis)
    associate(h=>c%sweep(d)%hydro(kk))
     allocate(h%volume(1-ng:n(d)+ng,1-ng:n(a)+ng),h%area(n(d)+1,n(a)), &
      h%frame(n(d)+1,n(a),9),h%background(n(d)+1,n(a),8))
     idx(baxis)=kk
     do jj=1-ng,n(a)+ng
      idx(a)=jj
      do ii=1-ng,n(d)+ng
       idx(d)=ii;h%volume(ii,jj)=c%vol(idx(1),idx(2),idx(3))
      enddo
     enddo
     do jj=1,n(a)
      idx(a)=jj
      do ii=1,n(d)+1
       idx(d)=ii;normal=c%area(:,idx(1),idx(2),idx(3),d);h%area(ii,jj)=norm2(normal)
       call frame_from_normal(normal,frame);h%frame(ii,jj,:)=frame
       ed=idx;ed(d)=ed(d)-1;h%background(ii,jj,:)=yy_hydro_full_quad(:,ed(1),ed(2),ed(3),d)
      enddo
     enddo
     c%sweep_cache_bytes=c%sweep_cache_bytes+8_int64*(size(h%volume,kind=int64)+size(h%area,kind=int64) &
      +size(h%frame,kind=int64)+size(h%background,kind=int64))
    end associate
   enddo
   do kk=1,n(d)
    associate(e=>c%sweep(d)%bfield(kk))
     allocate(e%xyz(1-ng:n(a)+ng,1-ng:n(baxis)+ng,3),e%area(1-ng:n(a)+ng,1-ng:n(baxis)+ng,2), &
      e%normal(1-ng:n(a)+ng,1-ng:n(baxis)+ng,3,2),e%projection(1-ng:n(a)+ng,1-ng:n(baxis)+ng,2), &
      e%length(n(a)+1,n(baxis)+1),e%frame(n(a)+1,n(baxis)+1,9),e%background(n(a)+1,n(baxis)+1,3))
     idx(d)=kk
     do jj=1-ng,n(baxis)+ng
      idx(baxis)=jj
      do ii=1-ng,n(a)+ng
       idx(a)=ii;ed=idx-1;high=ed;high(d)=high(d)+1
       e%xyz(ii,jj,:)=(c%xv(:,ed(1),ed(2),ed(3))+c%xv(:,high(1),high(2),high(3)))/2
       do cvar=1,2
        i=merge(a,baxis,cvar==1);normal=c%area(:,idx(1),idx(2),idx(3),i)
        e%area(ii,jj,cvar)=norm2(normal);e%normal(ii,jj,:,cvar)=normal/norm2(normal)
        ! Preserve the original dot-product/division expression exactly.
        e%projection(ii,jj,cvar)=dot_product(c%reference_uniform,normal)/norm2(normal)
       enddo
      enddo
     enddo
     do jj=1,n(baxis)+1
      idx(baxis)=jj
      do ii=1,n(a)+1
       idx(a)=ii;ed=idx-1;high=ed;high(d)=high(d)+1
       a0=c%xv(:,ed(1),ed(2),ed(3));b0=c%xv(:,high(1),high(2),high(3));edge=b0-a0
       e%length(ii,jj)=norm2(edge);frame(7:9)=edge/norm2(edge)
       low=ed;low(baxis)=low(baxis)-1;high=ed;high(baxis)=high(baxis)+1
       pm=c%xv(:,high(1),high(2),high(3))-c%xv(:,low(1),low(2),low(3))
       low(d)=low(d)+1;high(d)=high(d)+1
       pm=pm+c%xv(:,high(1),high(2),high(3))-c%xv(:,low(1),low(2),low(3))
       frame(1:3)=cube_cross(pm,frame(7:9));frame(1:3)=frame(1:3)/norm2(frame(1:3))
       frame(4:6)=cube_cross(frame(7:9),frame(1:3));e%frame(ii,jj,:)=frame
       ed=idx;ed(a)=ed(a)-1;ed(baxis)=ed(baxis)-1
       e%background(ii,jj,:)=c%background_edge(:,ed(1),ed(2),ed(3),d)+c%reference_uniform
      enddo
     enddo
     c%sweep_cache_bytes=c%sweep_cache_bytes+8_int64*(size(e%xyz,kind=int64)+size(e%area,kind=int64) &
      +size(e%normal,kind=int64)+size(e%projection,kind=int64)+size(e%length,kind=int64) &
      +size(e%frame,kind=int64)+size(e%background,kind=int64))
    end associate
   enddo
  enddo
  c%sweep_cache_ready=.true.
 end subroutine
 subroutine read_reference_background(path,rp,b)
  character(*),intent(in)::path
  real(real64),intent(in)::rp
  real(real64),allocatable,intent(out)::b(:,:,:,:)
  character(8)::magic
  integer(int32)::n(3),extra
  integer::u,st
  real(real64)::stored_rp,beq
  open(newunit=u,file=path,access='stream',form='unformatted',status='old',iostat=st)
  call demand(st==0,'MI static native background file')
  read(u,iostat=st) magic,n,stored_rp,beq
  call demand(st==0,'MI background header')
  call demand(magic=='LFMCB0N1'.and.all(n==[4,48,64]).and.stored_rp==rp.and.beq==.30_real64,'MI background convention/shape')
  allocate(b(4,48,64,3));read(u,iostat=st) b;call demand(st==0.and.all(ieee_is_finite(b)),'MI background finite payload')
  read(u,iostat=st) extra;call demand(st<0,'MI background trailing bytes');close(u)
 end subroutine
 subroutine model_couple(c)
  class(cube_model),intent(inout)::c
#ifdef LFM_CUBE_NATIVE_MI
  real(real64),allocatable::inner(:,:,:,:),planetary_phi(:,:,:,:)
  real(real64)::tick
  integer::st
  if(.not.c%mi_enabled.or..not.clock_reached(c%time,c%mi_next)) return
  tick=MPI_Wtime()
  ! MI samples the completed state, so synchronize its real current ghosts.
  call c%guards(c%q)
  call planetary_current_reference(c,planetary_phi)
  call yy_source_refresh(c%q,c%phi,planetary_phi,c%halo,ng,c%time,MPI_COMM_WORLD,st)
  call demand(st==0,'released perturbation curl about full planetary CT reference and IGRF projection')
  allocate(inner(4,c%n(2),c%n(3),8));inner=0
  if(c%layout%coord(1)==0) inner=c%q(1:4,1:c%n(2),1:c%n(3),:)
  ! Native packed state retains rho/Cs; FAC comes from the released CT source.
  call c%mi%update(inner,c%step,c%time,st,dt_s=max(0._real64,c%dt_old))
  call demand(st==0,'released Cube YY-current / IGRF mapping / LFM100 coupling')
  c%mi_updates=c%mi_updates+1
  do while(clock_reached(c%time,c%mi_next))
   c%mi_next=c%mi_next+c%mi_cadence
  enddo
  c%phase_seconds(7)=c%phase_seconds(7)+MPI_Wtime()-tick
#endif
 end subroutine
 subroutine planetary_current_reference(c,phi_ref)
  ! Same nativeGINT/curl face reference used by HYDRO; evaluate at the
  ! completed MI epoch, not the previous advance's cached reference time.
  ! The prescribed planetary field is curl-free in the MI source shell.
  ! Subtract its *discrete* CT representation before the unchanged donor
  ! recovery/curl. Actual perturbation flux and current are retained.
  class(cube_model),intent(inout)::c
  real(real64),allocatable,intent(out)::phi_ref(:,:,:,:)
  real(real64)::weights(21),phase
  integer::i,j,k,d,idx(3),ed(3),st
  call demand(allocated(yy_planetary_initial_phi),'planetary CT reference initialized')
  allocate(phi_ref,mold=c%phi)
  if(c%time==0._real64)then
   ! Reuse the exact known initial reference, avoiding coefficient-fit roundoff.
   phi_ref=yy_planetary_initial_phi
  else
   phi_ref=0;phase=yy_spin_omega_rad_s*c%time;weights(1)=1
   do i=1,10
    weights(2*i)=cos(i*phase);weights(2*i+1)=sin(i*phase)
   enddo
   do d=1,3
    do k=1,c%n(3)+merge(1,0,d==3)
     do j=1,c%n(2)+merge(1,0,d==2)
      do i=1,c%n(1)+merge(1,0,d==1)
       idx=[i,j,k];ed=idx;ed(d)=ed(d)-1
       phi_ref(ed(1),ed(2),ed(3),d)=dot_product(yy_hydro_quad_coef(1,:,ed(1),ed(2),ed(3),d),weights)*c%area_magnitude(i,j,k,d)
      enddo
     enddo
    enddo
   enddo
   call c%faces%sync(phi_ref,st);call demand(st==0,'MI planetary reference physical CT face ownership')
  endif
  call demand(all(ieee_is_finite(phi_ref)),'finite MI full planetary CT reference')
 end subroutine
 subroutine get_face(c,d,idx,v)
  class(cube_model),intent(in)::c
  integer,intent(in)::d,idx(3)
  real(real64),intent(out)::v(3,4)
  integer::a,b,q(3,4),z
  a=mod(d,3)+1;b=mod(d+1,3)+1;q=spread(idx-1,2,4);q(d,:)=idx(d)-1
  q(a,2:3)=idx(a);q(b,3:4)=idx(b)
  do z=1,4
   v(:,z)=c%xv(:,q(1,z),q(2,z),q(3,z))
  enddo
 end subroutine
 subroutine add_initial_dipole(c)
  class(cube_model),intent(inout)::c
  real(real64),allocatable::ap(:,:,:,:)
  real(real64)::a(3),b(3),s
  real(real64),external::gint,cube_yy_ax,cube_yy_ay,cube_yy_az
  integer::d,i,j,k,l,idx(3),p(3),st,n(3)
  call demand(cube_yy_background_ready(),'explicit tagged YY provider configuration required')
  call cube_yy_background_context(c%time,0)
  n=c%n;allocate(ap,mold=c%phi);ap=0
  ! Byte-extracted native GINT and cube_yy_ax/cube_yy_ay/cube_yy_az provide the original line rule.
  do d=1,3
   do k=0,n(3)
    do j=0,n(2)
     do i=0,n(1)
      idx=[i,j,k];if(idx(d)==0) cycle
      p=idx;p(d)=p(d)-1;a=c%xv(:,p(1),p(2),p(3));b=c%xv(:,i,j,k)
      s=(b(1)-a(1))*gint(cube_yy_ax,a(1),a(2),a(3),b(1),b(2),b(3)) &
       +(b(2)-a(2))*gint(cube_yy_ay,a(1),a(2),a(3),b(1),b(2),b(3)) &
       +(b(3)-a(3))*gint(cube_yy_az,a(1),a(2),a(3),b(1),b(2),b(3))
      ap(i,j,k,d)=s
     enddo
    enddo
   enddo
  enddo
  call c%edges%sync(ap,st);call demand(st==0,'initial shared vector potential')
  call cube_curl_update(c%background_flux,ap,n,-1._real64)
  c%phi=c%phi+c%background_flux
 end subroutine
 subroutine background_quadratures(c)
  class(cube_model),intent(inout)::c
  integer::i,j,k,d,idx(3),st,n(3),p(3)
  real(real64)::v(3,4),a(3),b(3),q(8)
  real(real64),external::gint,cube_yy_bxqq0,cube_yy_byqq0,cube_yy_bzqq0
  n=c%n
  do d=1,3
   do k=0,n(3)
    do j=0,n(2)
     do i=0,n(1)
      idx=[i,j,k]
      if(.not.any(pack(idx,[d/=1,d/=2,d/=3])==0)) then
       idx(d)=idx(d)+1;call get_face(c,d,idx,v)
       call cube_yy_g2int(q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),norm2(c%p%area(:,i,j,k,d)), &
        v(1,1),v(1,2),v(1,4),v(1,3),v(2,1),v(2,2),v(2,4),v(2,3), &
        v(3,1),v(3,2),v(3,4),v(3,3))
       c%background_quad(:,i,j,k,d)=q
       call cube_yy_background_pressure_delta(v,c%area(:,idx(1),idx(2),idx(3),d), &
        yy_background_pressure_correction(:,i,j,k,d))
       if(c%mode>=2) c%background_flux(i,j,k,d)=q(1)*norm2(c%p%area(:,i,j,k,d))
      endif
      idx=[i,j,k]
      if(idx(d)==0) cycle
      p=idx;p(d)=p(d)-1;a=c%xv(:,p(1),p(2),p(3));b=c%xv(:,i,j,k)
      c%background_edge(:,i,j,k,d)=[gint(cube_yy_bxqq0,a(1),a(2),a(3),b(1),b(2),b(3)), &
       gint(cube_yy_byqq0,a(1),a(2),a(3),b(1),b(2),b(3)),gint(cube_yy_bzqq0,a(1),a(2),a(3),b(1),b(2),b(3))]
     enddo
    enddo
   enddo
  enddo
  do k=1,n(3)
   do j=1,n(2)
    do i=1,n(1)
     c%background_cell(i,j,k,:)=(c%background_quad(3:5,i,j,k-1,3)+c%background_quad(3:5,i,j,k,3))/2
    enddo
   enddo
  enddo
 end subroutine
 ! Bind archived total B to the exact static background used by this case.
 ! Observational output only; no field calculation, clock or state change.
 subroutine prepare_hydro_time_reference(c)
  ! Same original CT representation as actual q: m0/G2 cell plus native BZZ
  ! of total flux minus stored m0/G2 normal flux. Do not substitute a
  ! different G2 face-average cell field for that native CT recovery.
  class(cube_model),intent(inout)::c
  integer,parameter::ns=21
  real(real64),allocatable::ap(:,:,:,:),phi_ref(:,:,:,:),delta_phi(:,:,:,:),recovered(:,:,:,:),reference_channels(:,:,:,:),reference_faces(:,:,:,:),edge_vector(:,:,:,:,:)
  real(real64)::v(3,4),qf(8),area,raw(11,ns),cells(3,ns),corr(3),a(3),b(3),line,t
  real(real64),external::gint,diagnostic_full_ax,diagnostic_full_ay,diagnostic_full_az,diagnostic_full_bx,diagnostic_full_by,diagnostic_full_bz
  integer::n(3),i,j,k,d,idx(3),ed(3),p(3),sample,st
  logical::bounded(8),positive(8)
  n=c%n
  call demand(.not.allocated(yy_hydro_full_quad),'one HYDRO reference per model')
  allocate(yy_hydro_full_quad,mold=c%background_quad);yy_hydro_full_quad=0
  allocate(yy_hydro_cell_shift(3,1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng));yy_hydro_cell_shift=0
  if(c%mode<1)return
  allocate(yy_hydro_quad_coef(11,ns,0:n(1),0:n(2),0:n(3),3));yy_hydro_quad_coef=0
  allocate(yy_hydro_cell_coef(3,ns,1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng))
  allocate(ap,mold=c%phi);allocate(phi_ref,mold=c%phi);allocate(delta_phi,mold=c%phi)
  allocate(recovered(n(1),n(2),n(3),3))
  allocate(reference_channels(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng,8))
  allocate(yy_bfield_face_coef(3,ns,1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng));yy_bfield_face_coef=0
  allocate(yy_bfield_face_reference(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng,3));yy_bfield_face_reference=0
  allocate(yy_bfield_edge_coef(3,ns,0:n(1),0:n(2),0:n(3),3));yy_bfield_edge_coef=0
  allocate(yy_bfield_edge_reference(3,0:n(1),0:n(2),0:n(3),3));yy_bfield_edge_reference=0
  allocate(edge_vector(3,0:n(1),0:n(2),0:n(3),3))
  bounded=.false.;positive=.false.
  do sample=1,ns
   t=yy_spin_period_s*(sample-1)/ns
   call diagnostic_full_set_time(t);call diagnostic_full_A_set_time(t)
   ap=0;phi_ref=0;edge_vector=0
   ! Exact original add_initial_dipole edge rule, GINT and MPI edge owner.
   do d=1,3
    do k=0,n(3)
     do j=0,n(2)
      do i=0,n(1)
       idx=[i,j,k];if(idx(d)==0)cycle
       p=idx;p(d)=p(d)-1;a=c%xv(:,p(1),p(2),p(3));b=c%xv(:,i,j,k)
       line=(b(1)-a(1))*gint(diagnostic_full_ax,a(1),a(2),a(3),b(1),b(2),b(3)) &
        +(b(2)-a(2))*gint(diagnostic_full_ay,a(1),a(2),a(3),b(1),b(2),b(3)) &
        +(b(3)-a(3))*gint(diagnostic_full_az,a(1),a(2),a(3),b(1),b(2),b(3))
       ap(i,j,k,d)=line
       edge_vector(:,i,j,k,d)=[gint(diagnostic_full_bx,a(1),a(2),a(3),b(1),b(2),b(3)), &
        gint(diagnostic_full_by,a(1),a(2),a(3),b(1),b(2),b(3)), &
        gint(diagnostic_full_bz,a(1),a(2),a(3),b(1),b(2),b(3))]
      enddo
     enddo
    enddo
   enddo
   call c%edges%sync(ap,st);call demand(st==0,'reference shared nativeGINT edge ownership')
   ! BG is a Cartesian vector average, as in original background_quadratures.
   ! Do not use the signed edge-integral owner map on these vector components.
   yy_bfield_edge_coef(:,sample,:,:,:,:)=edge_vector
   call cube_curl_update(phi_ref,ap,n,-1._real64)
   if(c%mode==3.and.sample==1)call demand(all(phi_ref==c%phi),'reference epoch0 actual CT flux equality')
   ! Apply the same physical face-owner reconciliation before BZZ recovery.
   call c%faces%sync(phi_ref,st);call demand(st==0,'reference physical CT face ownership')
   if(sample==1)then
    allocate(yy_planetary_initial_phi,mold=phi_ref);yy_planetary_initial_phi=phi_ref
   endif
   delta_phi=phi_ref-c%background_flux
   call native_recover(c%p,delta_phi,c%correction,recovered)
   ! Transport reference B through the SAME linear primitive halo. Other
   ! channels are unused packing slots, never a physical model state.
   reference_channels=0;reference_channels(1:n(1),1:n(2),1:n(3),6:8)=recovered
   call c%halo%sync(reference_channels,st,bounded,positive)
   call demand(st==0,'reference B native halo transport')
   do k=1-ng,n(3)+ng
    do j=1-ng,n(2)+ng
     if(c%layout%coord(1)==0)then
      do i=1-ng,0
       reference_channels(i,j,k,6:8)=reference_channels(1,j,k,6:8)
      enddo
     endif
     if(c%layout%coord(1)==c%layout%dims(1)-1)then
      ! Original prescribed planetary taper is identically zero here;
      ! nearest active planetary reference also gives the inflow value0.
      do i=n(1)+1,n(1)+ng
       reference_channels(i,j,k,6:8)=reference_channels(n(1),j,k,6:8)
      enddo
     endif
    enddo
   enddo
   call reference_staggered_fields(c,reference_channels,phi_ref,t,reference_faces)
   do k=1-ng,n(3)+ng
    do j=1-ng,n(2)+ng
     do i=1-ng,n(1)+ng
      yy_bfield_face_coef(:,sample,i,j,k)=reference_faces(i,j,k,:)
      yy_hydro_cell_coef(:,sample,i,j,k)=reference_channels(i,j,k,6:8)
     enddo
    enddo
   enddo
   deallocate(reference_faces)
   do d=1,3
    do k=1,n(3)+merge(1,0,d==3)
     do j=1,n(2)+merge(1,0,d==2)
      do i=1,n(1)+merge(1,0,d==1)
       idx=[i,j,k];ed=idx;ed(d)=ed(d)-1;call get_face(c,d,idx,v);area=norm2(cube_quad(v))
       call diagnostic_full_g2int(qf(1),qf(2),qf(3),qf(4),qf(5),qf(6),qf(7),qf(8),area, &
        v(1,1),v(1,2),v(1,4),v(1,3),v(2,1),v(2,2),v(2,4),v(2,3),v(3,1),v(3,2),v(3,4),v(3,3))
       yy_hydro_quad_coef(1:8,sample,ed(1),ed(2),ed(3),d)=qf
       yy_hydro_quad_coef(6:8,sample,ed(1),ed(2),ed(3),d)=qf(6:8)*(qf(1)+1.e-20_real64/area)
       yy_hydro_quad_coef(1,sample,ed(1),ed(2),ed(3),d)=phi_ref(ed(1),ed(2),ed(3),d)/area
       call diagnostic_full_pressure_delta(v,c%area(:,idx(1),idx(2),idx(3),d),corr)
       yy_hydro_quad_coef(9:11,sample,ed(1),ed(2),ed(3),d)=corr
      enddo
     enddo
    enddo
   enddo
   if(c%rank==0)write(*,*)'INIT_CT_MATCHED_HYDRO_REFERENCE_PHASE',sample,ns
  enddo
  deallocate(ap,phi_ref,delta_phi,recovered,reference_channels,edge_vector)
  do k=1-ng,n(3)+ng
   do j=1-ng,n(2)+ng
    do i=1-ng,n(1)+ng
     cells=yy_hydro_cell_coef(:,:,i,j,k)
     call yy_reference_coefficients(cells,yy_hydro_cell_coef(:,:,i,j,k))
     cells=yy_bfield_face_coef(:,:,i,j,k);call yy_reference_coefficients(cells,yy_bfield_face_coef(:,:,i,j,k))
    enddo
   enddo
  enddo
  do d=1,3
   do k=1,n(3)+merge(1,0,d==3)
    do j=1,n(2)+merge(1,0,d==2)
     do i=1,n(1)+merge(1,0,d==1)
      ed=[i,j,k];ed(d)=ed(d)-1;raw=yy_hydro_quad_coef(:,:,ed(1),ed(2),ed(3),d)
      call yy_reference_coefficients(raw,yy_hydro_quad_coef(:,:,ed(1),ed(2),ed(3),d))
     enddo
    enddo
   enddo
  enddo
  do d=1,3
   do k=0,n(3);do j=0,n(2);do i=0,n(1)
    cells=yy_bfield_edge_coef(:,:,i,j,k,d);call yy_reference_coefficients(cells,yy_bfield_edge_coef(:,:,i,j,k,d))
   enddo;enddo;enddo
  enddo
  call demand(all(ieee_is_finite(yy_bfield_face_coef)).and.all(ieee_is_finite(yy_bfield_edge_coef)),'finite fullplanet native BFIELD reference coefficients')
  call demand(all(ieee_is_finite(yy_hydro_quad_coef)).and.all(ieee_is_finite(yy_hydro_cell_coef)),'finite native CT matched reference coefficients')
  call diagnostic_full_set_time(0._real64);call diagnostic_full_A_set_time(0._real64)
  call prepare_zero_reference_records(c)
  call prepare_direct_reference_banks(c)
  call qualify_reference_evaluation(c)
  call update_hydro_time_reference(c)
  if(c%rank==0)write(*,*)'INIT_CT_MATCHED_HYDRO_REFERENCE_PASS; actual MHD CT/BFIELD/MI representation unchanged'
 end subroutine
 subroutine prepare_zero_reference_records(c)
  class(cube_model),intent(in)::c
  integer::n(3),i,j,k,d,ed(3)
  n=c%n
  allocate(yy_zero_cell(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng))
  allocate(yy_zero_face,mold=yy_zero_cell)
  allocate(yy_zero_edge(0:n(1),0:n(2),0:n(3),3))
  allocate(yy_zero_quad(0:n(1),0:n(2),0:n(3),3));yy_zero_quad=.true.
  do k=1-ng,n(3)+ng;do j=1-ng,n(2)+ng;do i=1-ng,n(1)+ng
   yy_zero_cell(i,j,k)=all(yy_hydro_cell_coef(:,:,i,j,k)==0._real64)
   yy_zero_face(i,j,k)=all(yy_bfield_face_coef(:,:,i,j,k)==0._real64)
  enddo;enddo;enddo
  do d=1,3
   do k=0,n(3);do j=0,n(2);do i=0,n(1)
    yy_zero_edge(i,j,k,d)=all(yy_bfield_edge_coef(:,:,i,j,k,d)==0._real64)
   enddo;enddo;enddo
   do k=1,n(3)+merge(1,0,d==3);do j=1,n(2)+merge(1,0,d==2);do i=1,n(1)+merge(1,0,d==1)
    ed=[i,j,k];ed(d)=ed(d)-1
    yy_zero_quad(ed(1),ed(2),ed(3),d)=all(yy_hydro_quad_coef(:,:,ed(1),ed(2),ed(3),d)==0._real64)
   enddo;enddo;enddo
  enddo
  write(*,*)'EXACT_ZERO_REFERENCE_COUNTS rank/cell/face/edge/quad=',c%rank, &
   count(yy_zero_cell),count(yy_zero_face),count(yy_zero_edge),count(yy_zero_quad)
 end subroutine
 ! Fixed extents avoid general MATMUL descriptor/allocation dispatch at
 ! every cell and face. Retain all21 terms, original FP64 coefficients,
 ! trig weights, and one ordered multiply/add per coefficient and output.
 pure subroutine reference_eval3(a,w,r)
  real(real64),intent(in)::a(3,21),w(21)
  real(real64),intent(out)::r(3)
  integer::j,k
  r=0
  do k=1,21
   do j=1,3
    r(j)=r(j)+a(j,k)*w(k)
   enddo
  enddo
 end subroutine
 pure subroutine reference_eval11(a,w,r)
  real(real64),intent(in)::a(11,21),w(21)
  real(real64),intent(out)::r(11)
  integer::j,k
  r=0
  do k=1,21
   do j=1,11
    r(j)=r(j)+a(j,k)*w(k)
   enddo
  enddo
 end subroutine
 subroutine direct_bank_new(bank,active,count)
  real(real64),allocatable,intent(out)::bank(:,:,:)
  integer(c_int),allocatable,intent(out)::active(:)
  integer,intent(in)::count
  allocate(bank(32,21,(count+31)/32),active((count+31)/32));bank=0;active=0
 end subroutine
 subroutine direct_bank_put(bank,index,a)
  real(real64),intent(inout)::bank(:,:,:)
  integer,intent(in)::index
  real(real64),intent(in)::a(:)
  call demand(size(a)==21,'same original21 coefficient column')
  bank(mod(index-1,32)+1,:,(index-1)/32+1)=a
 end subroutine
 subroutine direct_bank_finish(bank,active)
  real(real64),intent(in)::bank(:,:,:)
  integer(c_int),intent(out)::active(:)
  integer::b
  do b=1,size(active)
   active(b)=merge(1_c_int,0_c_int,any(bank(:,:,b)/=0._real64))
  enddo
 end subroutine
 subroutine prepare_direct_reference_banks(c)
  class(cube_model),intent(in)::c
  integer::n(3),nx,ny,nz,i,j,k,l,d,ed(3),ix,iy,iz,flat,ne
  n=c%n;nx=n(1)+2*ng;ny=n(2)+2*ng;nz=n(3)+2*ng
  call direct_bank_new(direct_cell_bank,direct_cell_active,3*nx*ny*nz)
  call direct_bank_new(direct_face_bank,direct_face_active,3*nx*ny*nz)
  ne=(n(1)+1)*(n(2)+1)*(n(3)+1)*3
  call direct_bank_new(direct_edge_bank,direct_edge_active,3*ne)
  call direct_bank_new(direct_quad_bank,direct_quad_active,11*ne)
  allocate(direct_quad_values(11,0:n(1),0:n(2),0:n(3),3));direct_quad_values=0
  do k=1-ng,n(3)+ng;do j=1-ng,n(2)+ng;do i=1-ng,n(1)+ng
   ix=i-(1-ng);iy=j-(1-ng);iz=k-(1-ng)
   do l=1,3
    flat=l+3*(ix+nx*(iy+ny*iz))
    call direct_bank_put(direct_cell_bank,flat,yy_hydro_cell_coef(l,:,i,j,k))
    flat=ix+nx*(iy+ny*(iz+nz*(l-1)))+1
    call direct_bank_put(direct_face_bank,flat,yy_bfield_face_coef(l,:,i,j,k))
   enddo
  enddo;enddo;enddo
  do d=1,3
   do k=0,n(3);do j=0,n(2);do i=0,n(1)
    do l=1,3
     flat=l+3*(i+(n(1)+1)*(j+(n(2)+1)*(k+(n(3)+1)*(d-1))))
     call direct_bank_put(direct_edge_bank,flat,yy_bfield_edge_coef(l,:,i,j,k,d))
    enddo
   enddo;enddo;enddo
   do k=1,n(3)+merge(1,0,d==3);do j=1,n(2)+merge(1,0,d==2);do i=1,n(1)+merge(1,0,d==1)
    ed=[i,j,k];ed(d)=ed(d)-1
    do l=1,11
     flat=l+11*(ed(1)+(n(1)+1)*(ed(2)+(n(2)+1)*(ed(3)+(n(3)+1)*(d-1))))
     call direct_bank_put(direct_quad_bank,flat,yy_hydro_quad_coef(l,:,ed(1),ed(2),ed(3),d))
    enddo
   enddo;enddo;enddo
  enddo
  call direct_bank_finish(direct_cell_bank,direct_cell_active)
  call direct_bank_finish(direct_face_bank,direct_face_active)
  call direct_bank_finish(direct_edge_bank,direct_edge_active)
  call direct_bank_finish(direct_quad_bank,direct_quad_active)
 end subroutine
 subroutine evaluate_direct_reference(w)
  real(real64),intent(in)::w(21)
  integer(c_int)::st
  st=direct_reference_eval(int(size(yy_hydro_cell_shift),c_int),direct_cell_bank,direct_cell_active,w,yy_hydro_cell_shift)
  call demand(st==0,'original21 sums direct CT cell layout')
  st=direct_reference_eval(int(size(yy_bfield_face_reference),c_int),direct_face_bank,direct_face_active,w,yy_bfield_face_reference)
  call demand(st==0,'original21 sums direct staggered face layout')
  st=direct_reference_eval(int(size(yy_bfield_edge_reference),c_int),direct_edge_bank,direct_edge_active,w,yy_bfield_edge_reference)
  call demand(st==0,'original21 sums direct edge layout')
  st=direct_reference_eval(int(size(direct_quad_values),c_int),direct_quad_bank,direct_quad_active,w,direct_quad_values)
  call demand(st==0,'original21 sums direct native quadrature layout')
 end subroutine
 subroutine reference_word_check(c,a,w,x,y,label)
  class(cube_model),intent(in)::c
  real(real64),intent(in)::a(:,:),w(:),x(:),y(:)
  character(*),intent(in)::label
  integer::l,u,ios
  logical::ok
  character(96)::name
  ok=.true.
  do l=1,size(x)
   if(transfer(x(l),0_int64)/=transfer(y(l),0_int64))ok=.false.
  enddo
  if(.not.ok)then
   write(name,'("reference-evaluation-mismatch-r",i3.3,".bin")')c%rank
   open(newunit=u,file=trim(name),access='stream',form='unformatted',status='new',iostat=ios)
   if(ios==0)then
    write(u)'YYREF013',int(shape(a),int32),int(size(w),int32),a,w,x,y
    close(u)
   endif
  endif
  call demand(ok,label)
 end subroutine
 subroutine qualify_reference_evaluation(c)
  class(cube_model),intent(in)::c
  real(real64),parameter::epochs(8)=[0._real64,.173456789_real64, &
   7200.0050561318985_real64,10800._real64,21600.104527889976_real64, &
   32400._real64,43199.823456789_real64,43200._real64]
  real(real64)::w(21),a3(3),b3(3),a11(11),b11(11),phase
  integer::n(3),i,j,k,d,l,s,ed(3)
  integer(int64)::checked
  if(.not.allocated(yy_hydro_quad_coef))return
  n=c%n;checked=0
  do s=1,size(epochs)
   phase=yy_spin_omega_rad_s*epochs(s);w(1)=1
   do i=1,10
    w(2*i)=cos(i*phase);w(2*i+1)=sin(i*phase)
   enddo
   call evaluate_direct_reference(w)
   do k=1-ng,n(3)+ng;do j=1-ng,n(2)+ng;do i=1-ng,n(1)+ng
    a3=matmul(yy_hydro_cell_coef(:,:,i,j,k),w)
    b3=yy_hydro_cell_shift(:,i,j,k)
    call reference_word_check(c,yy_hydro_cell_coef(:,:,i,j,k),w,a3,b3,'ordered3 CT cell reference WORD equality')
    a3=matmul(yy_bfield_face_coef(:,:,i,j,k),w)
    b3=yy_bfield_face_reference(i,j,k,:)
    call reference_word_check(c,yy_bfield_face_coef(:,:,i,j,k),w,a3,b3,'ordered3 staggered reference WORD equality')
    checked=checked+6
   enddo;enddo;enddo
   do d=1,3
    do k=0,n(3);do j=0,n(2);do i=0,n(1)
     a3=matmul(yy_bfield_edge_coef(:,:,i,j,k,d),w)
     b3=yy_bfield_edge_reference(:,i,j,k,d)
     call reference_word_check(c,yy_bfield_edge_coef(:,:,i,j,k,d),w,a3,b3,'ordered3 edge reference WORD equality')
     checked=checked+3
    enddo;enddo;enddo
    do k=1,n(3)+merge(1,0,d==3);do j=1,n(2)+merge(1,0,d==2);do i=1,n(1)+merge(1,0,d==1)
     ed=[i,j,k];ed(d)=ed(d)-1
     a11=matmul(yy_hydro_quad_coef(:,:,ed(1),ed(2),ed(3),d),w)
     b11=direct_quad_values(:,ed(1),ed(2),ed(3),d)
     call reference_word_check(c,yy_hydro_quad_coef(:,:,ed(1),ed(2),ed(3),d),w,a11,b11,'ordered11 HYDRO reference WORD equality')
     checked=checked+11
    enddo;enddo;enddo
   enddo
  enddo
  write(*,*)'EXACT_REFERENCE_EVALUATION_WORD_PASS rank/words=',c%rank,checked
 end subroutine
 subroutine update_hydro_time_reference(c)
  class(cube_model),intent(inout)::c
  real(real64)::values(11),shift(3),face_values(3),edge_values(3),phase,area,weights(21)
  integer::n(3),i,j,k,d,ed(3),idx(3),a,b,ii,jj,kk
  if(.not.allocated(yy_hydro_quad_coef))return
  n=c%n;phase=yy_spin_omega_rad_s*c%time
  ! Same qualified harmonic polynomial; trig weights computed once per time.
  weights(1)=1
  do i=1,10
   weights(2*i)=cos(i*phase);weights(2*i+1)=sin(i*phase)
  enddo
  call evaluate_direct_reference(weights)
  do d=1,3
   do k=1,n(3)+merge(1,0,d==3);do j=1,n(2)+merge(1,0,d==2);do i=1,n(1)+merge(1,0,d==1)
    idx=[i,j,k];ed=idx;ed(d)=ed(d)-1;area=c%area_magnitude(i,j,k,d)
    values=direct_quad_values(:,ed(1),ed(2),ed(3),d)
    yy_hydro_full_quad(:,ed(1),ed(2),ed(3),d)=values(1:8)
    yy_hydro_full_quad(6:8,ed(1),ed(2),ed(3),d)=values(6:8)/(values(1)+1.e-20_real64/area)
    yy_background_pressure_correction(:,ed(1),ed(2),ed(3),d)=values(9:11)
   enddo;enddo;enddo
  enddo
  if(c%sweep_cache_ready)then
   do d=1,3
    a=mod(d,3)+1;b=mod(d+1,3)+1
    do kk=1,n(b)
     idx(b)=kk
     associate(h=>c%sweep(d)%hydro(kk))
      do jj=1,n(a)
       idx(a)=jj
       do ii=1,n(d)+1
        idx(d)=ii;ed=idx;ed(d)=ed(d)-1
        h%background(ii,jj,:)=yy_hydro_full_quad(:,ed(1),ed(2),ed(3),d)
       enddo
      enddo
     end associate
    enddo
   enddo
  endif
 end subroutine
 subroutine export_static_background(c)
  class(cube_model),intent(in)::c
  integer::u,ios
  character(128)::name
  write(name,'("yy-static-background-r",i3.3,".bin")')c%rank
  open(newunit=u,file=trim(name),status='new',access='stream',form='unformatted',iostat=ios)
  call demand(ios==0,'new static-background provenance output')
  write(u,iostat=ios)'YYBGS001',int(c%rank,int32),int(c%n,int32),int(c%layout%face,int32), &
   int(c%layout%first,int32),c%rp,c%background_cell,c%background_flux
  call demand(ios==0,'complete already-computed static-background output')
  close(u,iostat=ios);call demand(ios==0,'closed static-background provenance output')
 end subroutine
 subroutine model_guards(c,q,predictor)
  class(cube_model),intent(inout)::c
  real(real64),intent(inout)::q(1-ng:,1-ng:,1-ng:,:)
  logical,optional,intent(in)::predictor
  logical::positive(8),bounded(8)
  real(real64)::normal(3),dotn,x(3),donor(8),tick,target_twice(3),ghost(3)
  integer::st,i,j,k,src,idx(3),n(3),v
  n=c%n;positive=.false.;positive([1,5])=.true.;bounded=.false.
  ! MHD_HACK extrapolates C1 without a sign restriction; HYDRO uses C1**2
  ! in its energy reconstruction, then clips against positive base-time C.
  ! C1 is a signed predictor, not a physical sound-speed state. Keep the
  ! density check and every base/updated/checkpoint sound-speed check.
  if(present(predictor)) then
   if(predictor) positive(5)=.false.
  endif
  call c%halo%sync(q,st,bounded,positive);call demand(st==0,'primitive halo')
  tick=MPI_Wtime()
  do k=1-ng,n(3)+ng
   do j=1-ng,n(2)+ng
    if(c%layout%coord(1)==0) then
     ! Reuse the released target once per physical wall face, not once per
     ! reflected ghost depth. The target is independent of source velocity.
     x=(c%xv(:,0,j-1,k-1)+c%xv(:,0,j,k-1)+c%xv(:,0,j,k)+c%xv(:,0,j-1,k))/4
     donor=q(1,j,k,:)
     call yy_boundary_ghost(x,c%time,donor(2:4),ghost,st)
     call demand(st==0,'released YY IGRF rotating gas wall')
     target_twice=ghost+donor(2:4)
     do i=1-ng,0
      src=1-i;q(i,j,k,:)=q(src,j,k,:)
      if(c%mode>=1) then
       q(i,j,k,2:4)=target_twice-q(src,j,k,2:4)
       ! Qualified donor GAMERA_STRICT_INNER_MAGNETIC_WALL: nearest B.
       q(i,j,k,6:8)=q(1,j,k,6:8)
      endif
     enddo
    endif
    if(c%layout%coord(1)==c%layout%dims(1)-1) then
     donor=q(n(1),j,k,:)
     do i=n(1)+1,n(1)+ng
      q(i,j,k,:)=q(n(1),j,k,:)
      if(c%mode>=1) then
       dotn=dot_product(c%sw(2:4),c%area(:,n(1)+1,j,k,1))
       if(dotn<0) q(i,j,k,:)=c%sw
      endif
      if(c%mode>=2) then
       x=c%outer_cell_xyz(:,i-n(1),j,k)
       normal=c%outer_normal(:,j,k)
       q(i,j,k,:)=kaiju_outer_state(c%sw,c%time,x,donor,normal,c%outer_cell_weight(i-n(1),j,k), &
        c%wind_switch_s,c%wind_bz_after,c%wind_reference_offset_s)
      endif
     enddo
    endif
   enddo
  enddo
  if(c%diagnostic_trace) then
   do v=1,5,4
    if(.not.positive(v).or.minval(q(:,:,:,v))>0) cycle
    idx=minloc(q(:,:,:,v))-ng
    write(*,*) 'GUARD_TRACE rank/step/component/local/global/time=',c%rank,c%step,v,idx, &
     c%layout%first+idx,c%time
    write(*,*) 'GUARD_TRACE value/history/current=',q(idx(1),idx(2),idx(3),v), &
     c%history(idx(1),idx(2),idx(3),v),c%q(idx(1),idx(2),idx(3),v)
   enddo
  endif
  call demand(all(ieee_is_finite(q)).and.minval(q(:,:,:,1))>0 &
   .and.(.not.positive(5).or.minval(q(:,:,:,5))>0),'finite positive guards')
  c%guard_boundary_seconds=c%guard_boundary_seconds+MPI_Wtime()-tick
 end subroutine
 subroutine cube_staggered_patch_guards(layout,bf,g,comm,status,profile)
  ! Exact same-face Cartesian-index copy, including the shared normal face.
  ! Pack both sides before either receive is installed. This retains the old
  ! pr=2 radial source payload and prevents same-axis updates feeding sends.
  ! Axis order1,2,3 propagates already-filled slabs to edge/corner guards.
  type(cube_layout),intent(in)::layout
  integer,intent(in)::g
  real(real64),intent(inout)::bf(1-g:,1-g:,1-g:,:)
  type(MPI_Comm),intent(in)::comm
  integer,intent(out)::status
  real(real64),optional,intent(inout)::profile(3)
  real(real64)::tick
  real(real64),allocatable::send_low(:,:,:,:),send_high(:,:,:,:),recv_low(:,:,:,:),recv_high(:,:,:,:)
  integer::n(3),axis,a,b,d,i,j,k,index(3),coord(3),low,high,st,amount,tag
  status=1;n=layout%local
  if(g<1.or.any(n<g+1).or.size(bf,4)/=3) return
  if(any([size(bf,1),size(bf,2),size(bf,3)]/=n+2*g)) return
  do axis=1,3
   if(layout%dims(axis)==1) cycle
   tick=MPI_Wtime()
   a=mod(axis,3)+1;b=mod(axis+1,3)+1
   low=MPI_PROC_NULL;high=MPI_PROC_NULL
   coord=layout%coord;coord(axis)=coord(axis)-1
   if(coord(axis)>=0) low=cube_rank(layout%face,coord,layout%dims)
   coord=layout%coord;coord(axis)=coord(axis)+1
   if(coord(axis)<layout%dims(axis)) high=cube_rank(layout%face,coord,layout%dims)
   allocate(send_low(g+1,1-g:n(a)+g,1-g:n(b)+g,3),send_high(g+1,1-g:n(a)+g,1-g:n(b)+g,3))
   allocate(recv_low(g+1,1-g:n(a)+g,1-g:n(b)+g,3),recv_high(g+1,1-g:n(a)+g,1-g:n(b)+g,3))
   do d=1,3
    do k=1-g,n(b)+g
     do j=1-g,n(a)+g
      index(a)=j;index(b)=k
      do i=1,g+1
       index(axis)=i;send_low(i,j,k,d)=bf(index(1),index(2),index(3),d)
       index(axis)=n(axis)-g+i;send_high(i,j,k,d)=bf(index(1),index(2),index(3),d)
      enddo
     enddo
    enddo
   enddo
   if(present(profile)) profile(1)=profile(1)+MPI_Wtime()-tick
   tick=MPI_Wtime()
   tag=171+2*(axis-1)
   call MPI_Sendrecv(send_high,size(send_high),MPI_DOUBLE_PRECISION,high,tag, &
    recv_low,size(recv_low),MPI_DOUBLE_PRECISION,low,tag,comm,MPI_STATUS_IGNORE,st)
   if(st/=MPI_SUCCESS) then
    status=st;return
   endif
   call MPI_Sendrecv(send_low,size(send_low),MPI_DOUBLE_PRECISION,low,tag+1, &
    recv_high,size(recv_high),MPI_DOUBLE_PRECISION,high,tag+1,comm,MPI_STATUS_IGNORE,st)
   if(st/=MPI_SUCCESS) then
    status=st;return
   endif
   if(present(profile)) profile(2)=profile(2)+MPI_Wtime()-tick
   tick=MPI_Wtime()
   do d=1,3
    do k=1-g,n(b)+g
     do j=1-g,n(a)+g
      index(a)=j;index(b)=k
      if(low/=MPI_PROC_NULL) then
       amount=g;if(d==axis) amount=g+1
       do i=1,amount
        index(axis)=i-g;bf(index(1),index(2),index(3),d)=recv_low(i,j,k,d)
       enddo
      endif
      if(high/=MPI_PROC_NULL) then
       do i=1,g
        index(axis)=n(axis)+i;bf(index(1),index(2),index(3),d)=recv_high(i,j,k,d)
       enddo
      endif
     enddo
    enddo
   enddo
   deallocate(send_low,send_high,recv_low,recv_high)
   if(present(profile)) profile(3)=profile(3)+MPI_Wtime()-tick
  enddo
  status=MPI_SUCCESS
 end subroutine
 subroutine model_face_guards(c,q,phi,bf,predictor_factor)
  class(cube_model),intent(inout)::c
  real(real64),intent(inout)::q(1-ng:,1-ng:,1-ng:,:)
  real(real64),intent(in)::phi(0:,0:,0:,:)
  real(real64),allocatable,intent(out)::bf(:,:,:,:)
  real(real64),optional,intent(in)::predictor_factor
  real(real64)::wall_tfac
  integer::i,j,k,d,idx(3),p(3),n(3),st,low(3),high(3)
  real(real64)::b(3),v(3,4),x(3),wind(8),tick,exact_wall_tick
  tick=MPI_Wtime();n=c%n;allocate(bf(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng,3))
  do d=1,3
   low=1;high=n;high(d)=high(d)+1
   do k=1-ng,n(3)+ng
    do j=1-ng,n(2)+ng
     do i=1-ng,n(1)+ng
      idx=[i,j,k]
      ! Active staggered faces are installed from the evolved CT integrals
      ! below. Their primitive reconstruction would be immediately overwritten.
      if(all(idx>=low).and.all(idx<=high)) cycle
      p=idx;p(d)=max(1-ng,p(d)-1)
      b=(q(i,j,k,6:8)+q(p(1),p(2),p(3),6:8))/2
      bf(i,j,k,d)=dot_product(b,c%area(:,i,j,k,d))/c%area_magnitude(i,j,k,d)
     enddo
    enddo
   enddo
  enddo
  c%face_guard_seconds(1)=c%face_guard_seconds(1)+MPI_Wtime()-tick;tick=MPI_Wtime()
  c%face_guard_seconds(2)=c%face_guard_seconds(2)+MPI_Wtime()-tick;tick=MPI_Wtime()
  do d=1,3
   ! Exactly the old valid index set and iteration order, expressed as loop
   ! bounds instead of an allocatable PACK result at every candidate face.
   low=1;low(d)=0
   do k=low(3),n(3)
    do j=low(2),n(2)
     do i=low(1),n(1)
      idx=[i,j,k];idx(d)=idx(d)+1
      bf(idx(1),idx(2),idx(3),d)=(phi(i,j,k,d)-c%background_flux(i,j,k,d))/c%area_magnitude(idx(1),idx(2),idx(3),d)
     enddo
    enddo
   enddo
  enddo
  c%face_guard_seconds(3)=c%face_guard_seconds(3)+MPI_Wtime()-tick
  tick=MPI_Wtime()
  call outer_face_fields(c,q,phi,bf)
  c%outer_seconds(1)=c%outer_seconds(1)+MPI_Wtime()-tick
  ! Same-patch normal/tangent face values are copied exactly on every tiled
  ! axis. Cross-patch angular guards retain the Cartesian reconstruction;
  ! evolving shared integrals remain owned by the existing CT plans.
  call cube_staggered_patch_guards(c%layout,bf,ng,MPI_COMM_WORLD,st,c%staggered_guard_seconds)
  call demand(st==MPI_SUCCESS,'exact same-patch staggered guards')
  if(c%mode>=1.and.c%layout%coord(1)==0) then
   wall_tfac=0;if(present(predictor_factor)) wall_tfac=predictor_factor
   exact_wall_tick=MPI_Wtime()
   call yy_magnetic_wall_apply(n(1),c%time,c%dt_old,wall_tfac,bf,st)
   equivalent_cost_s(5)=equivalent_cost_s(5)+MPI_Wtime()-exact_wall_tick
   call demand(st==0,'released rotating inner magnetic ghost flux')
  endif
 end subroutine
 ! A common numerical reconstruction reference from actual predicted state.
 ! Stored perturbation/CT fields and physical B0 remain unchanged. Reusing
 ! the native B0 interface avoids treating metric variation of a constant
 ! Cartesian field as a magnetic gradient. One reference across every tile
 ! avoids a tile-dependent field split. It is derived, not restart state.
 subroutine reference_staggered_fields(c,q,phi,time,bf)
  class(cube_model),intent(inout)::c
  real(real64),intent(in)::q(1-ng:,1-ng:,1-ng:,:),phi(0:,0:,0:,:),time
  real(real64),allocatable,intent(out)::bf(:,:,:,:)
  integer::n(3),i,j,k,d,idx(3),p(3),low(3),high(3),st
  real(real64)::b(3)
  n=c%n;allocate(bf(1-ng:n(1)+ng,1-ng:n(2)+ng,1-ng:n(3)+ng,3))
  do d=1,3
   low=1;high=n;high(d)=high(d)+1
   do k=1-ng,n(3)+ng;do j=1-ng,n(2)+ng;do i=1-ng,n(1)+ng
    idx=[i,j,k];if(all(idx>=low).and.all(idx<=high))cycle
    p=idx;p(d)=max(1-ng,p(d)-1);b=(q(i,j,k,6:8)+q(p(1),p(2),p(3),6:8))/2
    bf(i,j,k,d)=dot_product(b,c%area(:,i,j,k,d))/c%area_magnitude(i,j,k,d)
   enddo;enddo;enddo
   low=1;low(d)=0
   do k=low(3),n(3);do j=low(2),n(2);do i=low(1),n(1)
    idx=[i,j,k];idx(d)=idx(d)+1
    bf(idx(1),idx(2),idx(3),d)=(phi(i,j,k,d)-c%background_flux(i,j,k,d))/c%area_magnitude(idx(1),idx(2),idx(3),d)
   enddo;enddo;enddo
  enddo
  ! Prescribed planetary reference is identically0 at the physical outer
  ! boundary under the existing16..24 taper; solar wind is a perturbation.
  call cube_staggered_patch_guards(c%layout,bf,ng,MPI_COMM_WORLD,st)
  call demand(st==MPI_SUCCESS,'reference original same-patch staggered guards')
  if(c%mode>=1.and.c%layout%coord(1)==0)then
   call yy_magnetic_wall_apply(n(1),time,0._real64,0._real64,bf,st)
   call demand(st==0,'reference released rotating magnetic wall operator')
  endif
 end subroutine
 subroutine cube_magnetic_reference(c,reference)
  class(cube_model),intent(in)::c
  real(real64),intent(out)::reference(3)
  real(real64)::local(4),global(4),weight
  integer::i,j,k,st
  local=0
  do k=1,c%n(3)
   do j=1,c%n(2)
    do i=1,c%n(1)
     weight=c%vol(i,j,k)
     local(1:3)=local(1:3)+weight*(c%pred(i,j,k,6:8)-yy_hydro_cell_shift(:,i,j,k))
     local(4)=local(4)+weight
    enddo
   enddo
  enddo
  call MPI_Allreduce(local,global,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,st)
  call demand(st==MPI_SUCCESS.and.global(4)>0.and.all(ieee_is_finite(global)),'finite common magnetic reference')
  reference=global(1:3)/global(4)
 end subroutine
 subroutine seam_base_limiter_reference(c,d,kk,qb)
  class(cube_model),intent(in)::c
  integer,intent(in)::d,kk
  real(real64),intent(inout)::qb(1-ng:,1-ng:,:)
  integer::a,baxis,j,side,idx(3),h,ghost
  a=mod(d,3)+1;baxis=mod(d+1,3)+1
  do side=0,1
   idx(d)=side*c%n(d);ghost=merge(0,c%n(d)+1,side==0)
   idx(baxis)=kk
   do j=1,c%n(a)
    idx(a)=j
    h=c%fluxes%cell_pair_lookup(idx(1),idx(2),idx(3),d)
    call demand(h>0,'complete physical face BASE-reference lookup')
    if(.not.c%fluxes%cell_pair_mask(h)) cycle
    qb(ghost,j,1:5)=c%hydro_base_reference(1:5,h)
   enddo
  enddo
 end subroutine
 subroutine plane_word_check(c,q,p,rq,rp,label,d,k)
  class(cube_model),intent(in)::c
  real(real64),intent(in)::q(:,:,:),p(:,:,:),rq(:,:,:),rp(:,:,:)
  character(*),intent(in)::label
  integer,intent(in)::d,k
  integer::i,j,l,u,ios
  logical::ok
  character(160)::name
  call demand(all(shape(q)==shape(rq)).and.all(shape(p)==shape(rp)),'same actual plane oracle shape')
  ok=.true.
  do l=1,size(q,3);do j=1,size(q,2);do i=1,size(q,1)
   if(transfer(q(i,j,l),0_int64)/=transfer(rq(i,j,l),0_int64))ok=.false.
   if(transfer(p(i,j,l),0_int64)/=transfer(rp(i,j,l),0_int64))ok=.false.
  enddo;enddo;enddo
  if(.not.ok)then
   write(name,'("plane-word-mismatch-",a,"-r",i3.3,"-s",i7.7,"-d",i1,"-k",i3.3,".bin")')label,c%rank,c%step,d,k
   open(newunit=u,file=trim(name),status='new',access='stream',form='unformatted',iostat=ios)
   if(ios==0)then
    write(u)'YYPACK36',int(shape(q),int32),int([c%rank,c%step,d,k],int32),c%time,q,rq,p,rp
    close(u)
   endif
  endif
  call demand(ok,'actual original and contiguous plane WORD equality')
  plane_words_compared=plane_words_compared+2_int64*size(q,kind=int64)
 end subroutine
 subroutine model_advance(c,stop,dt,errors)
  class(cube_model),intent(inout)::c
  real(real64),intent(in)::stop
  real(real64),intent(out)::dt,errors(3)
  real(real64),allocatable::f(:,:,:,:,:),emf(:,:,:,:),bf(:,:,:,:),bp(:,:,:,:),phip(:,:,:,:),change(:,:,:,:)
  real(real64),allocatable::q0(:,:,:,:),before(:,:,:),audit(:,:,:,:)
  real(real64),allocatable::qb(:,:,:),qp(:,:,:),vv(:,:),aa(:,:),tt(:,:,:),bn(:,:),bg(:,:,:),ff(:,:,:),ac(:,:,:)
  real(real64),allocatable::xyz(:,:,:),length(:,:),frames(:,:,:),areas(:,:,:),normals(:,:,:,:),bb(:,:,:),pp(:,:,:)
  real(real64),allocatable::background(:,:,:),ee(:,:)
  real(real64),allocatable::copy_qb(:,:,:),copy_qp(:,:,:),copy_bb(:,:,:),copy_pp(:,:,:)
  logical::copy_fast,copy_check
#ifdef LFM_CUBE_NATIVE_MI
  real(real64),allocatable::wall_phi(:,:)
#endif
  real(real64)::tick_boundary,rotation_emf,edge_a(3),edge_b(3),exact_cost_tick,exact_emf_tick
  real(real64)::emf_reference(3),projected_reference
  real(real64)::localdt,tfac,frame(9),v4(3,4),normal(3),edge(3),a0(3),b0(3),p0(3),pm(3),mass,eng,u(3),rho
  real(real64)::mom(3),magmom(3),b(3),bdot,rb,sound2,scale,jump,wholejump,denom,wind(8),mix,tick,advance_start
  integer::n(3),d,a,baxis,i,j,k,cvar,ii,jj,kk,idx(3),ed(3),low(3),high(3),st,failure_unit,failure_io,rotation_index
  character(96)::failure_name
  n=c%n;advance_start=MPI_Wtime();tick=advance_start
  copy_fast=c%sweep_cache_enabled.and.c%mode>=2.and.all(n==[16,12,12])
  copy_check=c%step<8.or.iand(c%step,c%step-1)==0
  exact_cost_tick=MPI_Wtime()
  if(rotation_work_ready)then
   st=distributed_wall_prefetch(int(rotation_work_comm%MPI_VAL,c_int),c%time,c%dt_old)
   call demand(st==0,'original magnetic-wall current and original previous epoch workers')
  endif
  call refresh_rotation_work(c%time)
  equivalent_cost_s(8)=equivalent_cost_s(8)+MPI_Wtime()-exact_cost_tick
  call demand(.not.c%sweep_cache_enabled.or.c%sweep_cache_ready,'immutable sweep cache initialized')
  exact_cost_tick=MPI_Wtime()
  call c%guards(c%q)
  equivalent_cost_s(1)=equivalent_cost_s(1)+MPI_Wtime()-exact_cost_tick
  c%phase_seconds(1)=c%phase_seconds(1)+MPI_Wtime()-tick;tick=MPI_Wtime()
  localdt=native_cfl(c%q(1:n(1),1:n(2),1:n(3),:),n,ca,c%background_quad(2,:,:,:,3))
  call MPI_Allreduce(localdt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,st)
  c%cfl_candidate_dt=.3_real64*dt
  dt=min(c%cfl_candidate_dt,6._real64,stop-c%time);call demand(dt>0.and.ieee_is_finite(dt),'positive CFL step')
  c%phase_seconds(2)=c%phase_seconds(2)+MPI_Wtime()-tick;tick=MPI_Wtime()
  tfac=0;if(c%dt_old>0) tfac=.5_real64*dt/c%dt_old
  q0=c%q;c%pred=(1+tfac)*c%q-tfac*c%history
  exact_cost_tick=MPI_Wtime()
  call c%guards(c%pred,predictor=.true.)
  equivalent_cost_s(2)=equivalent_cost_s(2)+MPI_Wtime()-exact_cost_tick
  phip=(1+tfac)*c%phi-tfac*c%phi_history
  exact_cost_tick=MPI_Wtime()
  call c%face_guards(c%q,c%phi,bf);call c%face_guards(c%pred,phip,bp,predictor_factor=tfac)
  equivalent_cost_s(3)=equivalent_cost_s(3)+MPI_Wtime()-exact_cost_tick
  exact_cost_tick=MPI_Wtime()
  call update_hydro_time_reference(c)
  equivalent_cost_s(4)=equivalent_cost_s(4)+MPI_Wtime()-exact_cost_tick
  c%phase_seconds(1)=c%phase_seconds(1)+MPI_Wtime()-tick;tick=MPI_Wtime()
  call c%fluxes%exchange_base(ng,c%q,c%hydro_base_reference,st)
  call demand(st==0,'finite positive actual adjacent BASE references')
  allocate(f(0:n(1),0:n(2),0:n(3),3,8),emf(0:n(1),0:n(2),0:n(3),3));f=0;emf=0
  do d=1,3
   a=mod(d,3)+1;baxis=mod(d+1,3)+1
   allocate(qb(1-ng:n(d)+ng,1-ng:n(a)+ng,8),qp(1-ng:n(d)+ng,1-ng:n(a)+ng,8), &
    bn(n(d)+1,n(a)),ff(n(d)+1,n(a),8),ac(n(d),n(a),8))
   if(.not.c%sweep_cache_enabled) then
    allocate(vv(1-ng:n(d)+ng,1-ng:n(a)+ng),aa(n(d)+1,n(a)),tt(n(d)+1,n(a),9),bg(n(d)+1,n(a),8));bg=0
   endif
   do kk=1,n(baxis)
    idx(baxis)=kk
    if(.not.copy_fast.or.copy_check)then
    do jj=1-ng,n(a)+ng
     idx(a)=jj
     do ii=1-ng,n(d)+ng
      idx(d)=ii;qb(ii,jj,:)=c%q(idx(1),idx(2),idx(3),:);qp(ii,jj,:)=c%pred(idx(1),idx(2),idx(3),:)
      qb(ii,jj,6:8)=qb(ii,jj,6:8)-yy_hydro_cell_shift(:,idx(1),idx(2),idx(3))
      qp(ii,jj,6:8)=qp(ii,jj,6:8)-yy_hydro_cell_shift(:,idx(1),idx(2),idx(3))
      if(.not.c%sweep_cache_enabled) vv(ii,jj)=c%vol(idx(1),idx(2),idx(3))
     enddo
    enddo
     if(copy_fast)then
      copy_qb=qb;copy_qp=qp
     endif
    endif
    if(copy_fast)then
     st=pack_hydro_plane(n,ng,d,kk,c%q,c%pred,yy_hydro_cell_shift,qb,qp)
     call demand(st==0,'qualified contiguous HYDRO plane packing')
     if(copy_check)call plane_word_check(c,qb,qp,copy_qb,copy_qp,'HYDRO',d,kk)
    endif
    call seam_base_limiter_reference(c,d,kk,qb)
    do jj=1,n(a)
     idx(a)=jj
     do ii=1,n(d)+1
      idx(d)=ii;ed=idx;ed(d)=ed(d)-1
      bn(ii,jj)=bp(idx(1),idx(2),idx(3),d) &
       -(yy_hydro_full_quad(1,ed(1),ed(2),ed(3),d)-c%background_quad(1,ed(1),ed(2),ed(3),d))
      if(.not.c%sweep_cache_enabled) then
       normal=c%area(:,idx(1),idx(2),idx(3),d);aa(ii,jj)=norm2(normal)
       call frame_from_normal(normal,frame);tt(ii,jj,:)=frame
       ed=idx;ed(d)=ed(d)-1;bg(ii,jj,:)=yy_hydro_full_quad(:,ed(1),ed(2),ed(3),d)
      endif
     enddo
    enddo
    if(c%sweep_cache_enabled) then
     associate(h=>c%sweep(d)%hydro(kk))
      call lfm_cube_hydro_flux(n(d),n(a),ng,0._real64,gamma,ca,.75_real64,qb,qp, &
       h%volume,h%area,h%frame,bn,h%background,ff,ac,st)
     end associate
    else
     call lfm_cube_hydro_flux(n(d),n(a),ng,0._real64,gamma,ca,.75_real64,qb,qp,vv,aa,tt,bn,bg,ff,ac,st)
    endif
    call demand(st==0,'native HYDRO directional sweep')
    do jj=1,n(a)
     idx(a)=jj
     do ii=1,n(d)+1
      idx(d)=ii-1;f(idx(1),idx(2),idx(3),d,:)=ff(ii,jj,:)
      f(idx(1),idx(2),idx(3),d,6:8)=f(idx(1),idx(2),idx(3),d,6:8)+ &
       yy_background_pressure_correction(:,idx(1),idx(2),idx(3),d)
     enddo
    enddo
   enddo
   deallocate(qb,qp,bn,ff,ac)
   if(.not.c%sweep_cache_enabled) deallocate(vv,aa,tt,bg)
  enddo
  ! Released YY physical inner fluid-flux hook. Same +radial sign as the
  ! native face vector: suppress injection only; momentum/Maxwell unchanged.
  if(c%mode>=1.and.c%layout%coord(1)==0) then
   st=yy_flux_trap(n(1),n(2),n(3),f)
   call demand(st==0,'released inner outward mass energy flux boundary')
  endif
  call c%fluxes%sync_fields(f,st);call demand(st==0,'shared hydro face fluxes')
  allocate(change(n(1),n(2),n(3),8))
  do cvar=1,8
   change(:,:,:,cvar)=-cube_divergence(f(:,:,:,:,cvar),n)
  enddo
  before=cube_divergence(c%phi,n)
  c%phase_seconds(3)=c%phase_seconds(3)+MPI_Wtime()-tick;tick=MPI_Wtime()
  exact_emf_tick=MPI_Wtime()
  emf_reference=0
  if(c%mode>=2) call cube_magnetic_reference(c,emf_reference)
  do d=1,3
   a=mod(d,3)+1;baxis=mod(d+1,3)+1
   allocate(qb(1-ng:n(a)+ng,1-ng:n(baxis)+ng,8),qp(1-ng:n(a)+ng,1-ng:n(baxis)+ng,8), &
    bb(1-ng:n(a)+ng,1-ng:n(baxis)+ng,2),pp(1-ng:n(a)+ng,1-ng:n(baxis)+ng,2),ee(n(a)+1,n(baxis)+1))
   if(c%sweep_cache_enabled.and.c%mode>=2) allocate(background(n(a)+1,n(baxis)+1,3))
   if(.not.c%sweep_cache_enabled) then
    allocate(xyz(1-ng:n(a)+ng,1-ng:n(baxis)+ng,3),areas(1-ng:n(a)+ng,1-ng:n(baxis)+ng,2), &
     normals(1-ng:n(a)+ng,1-ng:n(baxis)+ng,3,2),background(n(a)+1,n(baxis)+1,3), &
     length(n(a)+1,n(baxis)+1),frames(n(a)+1,n(baxis)+1,9))
    do cvar=1,3
     background(:,:,cvar)=c%reference_uniform(cvar)
    enddo
   endif
   do kk=1,n(d)
    idx(d)=kk
    if(.not.copy_fast.or.copy_check)then
    do jj=1-ng,n(baxis)+ng
     idx(baxis)=jj
     do ii=1-ng,n(a)+ng
      idx(a)=ii;qb(ii,jj,:)=c%q(idx(1),idx(2),idx(3),:);qp(ii,jj,:)=c%pred(idx(1),idx(2),idx(3),:)
      if(c%sweep_cache_enabled) then
       do cvar=1,2
        i=merge(a,baxis,cvar==1)
        projected_reference=c%sweep(d)%bfield(kk)%projection(ii,jj,cvar)
        if(c%mode>=2) projected_reference=yy_bfield_face_reference(idx(1),idx(2),idx(3),i) &
         +dot_product(c%reference_uniform+emf_reference,c%area(:,idx(1),idx(2),idx(3),i))/c%area_magnitude(idx(1),idx(2),idx(3),i)
        bb(ii,jj,cvar)=bf(idx(1),idx(2),idx(3),i)-projected_reference
        pp(ii,jj,cvar)=bp(idx(1),idx(2),idx(3),i)-projected_reference
       enddo
      else
       ed=idx-1;high=ed;high(d)=high(d)+1
       xyz(ii,jj,:)=(c%xv(:,ed(1),ed(2),ed(3))+c%xv(:,high(1),high(2),high(3)))/2
       do cvar=1,2
        i=merge(a,baxis,cvar==1);normal=c%area(:,idx(1),idx(2),idx(3),i)
        areas(ii,jj,cvar)=norm2(normal);normals(ii,jj,:,cvar)=normal/norm2(normal)
        ! Reuse native BFIELD static-background representation for uniform
        ! IMF. PDM reconstructs perturbations; exact Cartesian background is
        ! supplied at edges, preserving uniform B on a curvilinear mesh.
        bb(ii,jj,cvar)=bf(idx(1),idx(2),idx(3),i)-dot_product(c%reference_uniform,normal)/norm2(normal)
        if(c%mode>=2) bb(ii,jj,cvar)=bf(idx(1),idx(2),idx(3),i) &
         -yy_bfield_face_reference(idx(1),idx(2),idx(3),i)-dot_product(c%reference_uniform+emf_reference,normal)/norm2(normal)
        pp(ii,jj,cvar)=bp(idx(1),idx(2),idx(3),i)-dot_product(c%reference_uniform,normal)/norm2(normal)
        if(c%mode>=2) pp(ii,jj,cvar)=bp(idx(1),idx(2),idx(3),i) &
         -yy_bfield_face_reference(idx(1),idx(2),idx(3),i)-dot_product(c%reference_uniform+emf_reference,normal)/norm2(normal)
       enddo
      endif
     enddo
    enddo
     if(copy_fast)then
      copy_qb=qb;copy_qp=qp;copy_bb=bb;copy_pp=pp
     endif
    endif
    if(copy_fast)then
     st=pack_bfield_plane(n,ng,d,kk,c%q,c%pred,bf,bp,yy_bfield_face_reference,c%area,c%area_magnitude, &
      c%reference_uniform,emf_reference,qb,qp,bb,pp)
     call demand(st==0,'qualified contiguous BFIELD plane packing')
     if(copy_check)then
      call plane_word_check(c,qb,qp,copy_qb,copy_qp,'BFIELD_Q',d,kk)
      call plane_word_check(c,bb,pp,copy_bb,copy_pp,'BFIELD_FACE',d,kk)
     endif
    endif
    if(.not.c%sweep_cache_enabled) then
     do jj=1,n(baxis)+1
      idx(baxis)=jj
      do ii=1,n(a)+1
       idx(a)=ii;ed=idx-1;high=ed;high(d)=high(d)+1
       a0=c%xv(:,ed(1),ed(2),ed(3));b0=c%xv(:,high(1),high(2),high(3));edge=b0-a0
       length(ii,jj)=norm2(edge);frame(7:9)=edge/norm2(edge)
       ! Original metric's edge-tangent plane, evaluated on actual vertices.
       low=ed;low(baxis)=low(baxis)-1;high=ed;high(baxis)=high(baxis)+1
       pm=c%xv(:,high(1),high(2),high(3))-c%xv(:,low(1),low(2),low(3))
       low(d)=low(d)+1;high(d)=high(d)+1
       pm=pm+c%xv(:,high(1),high(2),high(3))-c%xv(:,low(1),low(2),low(3))
       frame(1:3)=cube_cross(pm,frame(7:9));frame(1:3)=frame(1:3)/norm2(frame(1:3))
       frame(4:6)=cube_cross(frame(7:9),frame(1:3));frames(ii,jj,:)=frame
       ed=idx;ed(a)=ed(a)-1;ed(baxis)=ed(baxis)-1
       background(ii,jj,:)=c%background_edge(:,ed(1),ed(2),ed(3),d)+c%reference_uniform
       if(c%mode>=2) background(ii,jj,:)=yy_bfield_edge_reference(:,ed(1),ed(2),ed(3),d)+c%reference_uniform+emf_reference
      enddo
     enddo
    endif
    if(c%sweep_cache_enabled) then
     associate(e=>c%sweep(d)%bfield(kk))
      if(c%mode>=2) then
       do jj=1,n(baxis)+1
        idx(baxis)=jj-1
        do ii=1,n(a)+1
         idx(a)=ii-1
         background(ii,jj,:)=yy_bfield_edge_reference(:,idx(1),idx(2),idx(3),d)+c%reference_uniform+emf_reference
        enddo
       enddo
       call lfm_cube_bfield(n(a),n(baxis),ng,qb,qp,e%xyz,e%length,e%frame,e%area,e%normal,bb,pp,background,ee,st)
      else
      call lfm_cube_bfield(n(a),n(baxis),ng,qb,qp,e%xyz,e%length,e%frame,e%area,e%normal,bb,pp,e%background,ee,st)
      endif
     end associate
    else
     call lfm_cube_bfield(n(a),n(baxis),ng,qb,qp,xyz,length,frames,areas,normals,bb,pp,background,ee,st)
    endif
    call demand(st==0,'native BFIELD edge sweep')
    do jj=1,n(baxis)+1
     idx(baxis)=jj-1
     do ii=1,n(a)+1
      idx(a)=ii-1;emf(idx(1),idx(2),idx(3),d)=ee(ii,jj)
     enddo
    enddo
   enddo
   deallocate(qb,qp,bb,pp,ee)
   if(.not.c%sweep_cache_enabled) deallocate(xyz,areas,normals,background,length,frames)
   if(c%sweep_cache_enabled.and.c%mode>=2) deallocate(background)
  enddo
  equivalent_cost_s(6)=equivalent_cost_s(6)+MPI_Wtime()-exact_emf_tick
  tick_boundary=MPI_Wtime()
  call kaiju_outer_emf(c,c%q,dt,emf)
  c%outer_seconds(2)=c%outer_seconds(2)+MPI_Wtime()-tick_boundary
  exact_emf_tick=MPI_Wtime()
  if((c%mode>=1).and.c%layout%coord(1)==0) then
   emf(0,:,:,2:3)=0
#ifdef LFM_CUBE_NATIVE_MI
   if(c%mi_enabled) then
    call c%mi%wall_values(wall_phi,st,c%time);call demand(st==0,'held ionosphere wall voltage')
    rotation_index=0
    do k=0,n(3)
     do j=1,n(2)
      edge_a=.01_real64*c%xv(:,0,j-1,k);edge_b=.01_real64*c%xv(:,0,j,k)
      rotation_index=rotation_index+1
      if(rotation_work_ready.and.rotation_work_time==c%time)then
       rotation_emf=rotation_edge_values(rotation_index);st=0
      else
       call yy_corotation_edge(edge_a,edge_b,c%time,rotation_emf,st)
      endif
      call demand(st==0,'released YY rotating inner J edge')
      emf(0,j,k,2)=mi_edge_lfm(wall_phi(j-1,k),wall_phi(j,k))+rotation_emf
     enddo
    enddo
    do k=1,n(3)
     do j=0,n(2)
      edge_a=.01_real64*c%xv(:,0,j,k-1);edge_b=.01_real64*c%xv(:,0,j,k)
      rotation_index=rotation_index+1
      if(rotation_work_ready.and.rotation_work_time==c%time)then
       rotation_emf=rotation_edge_values(rotation_index);st=0
      else
       call yy_corotation_edge(edge_a,edge_b,c%time,rotation_emf,st)
      endif
      call demand(st==0,'released YY rotating inner K edge')
      emf(0,j,k,3)=mi_edge_lfm(wall_phi(j,k-1),wall_phi(j,k))+rotation_emf
     enddo
    enddo
   endif
#endif
  endif
  equivalent_cost_s(7)=equivalent_cost_s(7)+MPI_Wtime()-exact_emf_tick
  call c%edges%sync(emf,st);call demand(st==0,'shared native edge EMFs')
  c%phase_seconds(4)=c%phase_seconds(4)+MPI_Wtime()-tick;tick=MPI_Wtime()
  c%phi_history=c%phi
  call native_ct(c%phi,emf,n,dt)
  ! Original MHD_HACK accumulator -> primitive/Boris algebra, same order and
  ! native sound-speed floor1 cm/s. No density floor or velocity clamp added.
  do k=1,n(3)
   do j=1,n(2)
    do i=1,n(1)
     mass=q0(i,j,k,1)*c%vol(i,j,k)+dt*change(i,j,k,1)
     call demand(mass>0,'positive updated mass')
     eng=q0(i,j,k,1)*c%vol(i,j,k)*(sum(q0(i,j,k,2:4)**2)/2+q0(i,j,k,5)**2/(gamma*(gamma-1))) &
      +dt*change(i,j,k,5)
     mom=dt*change(i,j,k,2:4);magmom=dt*change(i,j,k,6:8)
     u=(mom+q0(i,j,k,2:4)*q0(i,j,k,1)*c%vol(i,j,k))/mass
     sound2=gamma*(gamma-1)*(eng/mass-sum(u*u)/2)
     if(c%diagnostic_trace.and.sound2<=1) then
      write(*,*) 'SOUND_TRACE rank/step/local/global/time=',c%rank,c%step,[i,j,k], &
       c%layout%first+[i,j,k],c%time
      write(*,*) 'SOUND_TRACE raw/eng/mass/u/q0=',sound2,eng,mass,u,q0(i,j,k,1:5)
     endif
     sound2=max(1._real64,sound2)
     c%q(i,j,k,5)=sqrt(max(1._real64,sound2));rho=mass/c%vol(i,j,k);c%q(i,j,k,1)=rho
     mom=mom-c%vol(i,j,k)*(rho-q0(i,j,k,1))*q0(i,j,k,2:4)
     b=c%pred(i,j,k,6:8)+c%background_cell(i,j,k,:);bdot=dot_product(b,mom)/(4*acos(-1._real64)*ca**2*rho)
     rb=c%vol(i,j,k)*(rho+sum(b*b)/(4*acos(-1._real64)*ca**2))
     c%q(i,j,k,2:4)=q0(i,j,k,2:4)+(mom+magmom+b*bdot)/rb
    enddo
   enddo
  enddo
  call native_recover(c%p,c%phi-c%background_flux,c%correction,c%q(1:n(1),1:n(2),1:n(3),6:8))
  tick_boundary=MPI_Wtime()
  call kaiju_nudge_outer(c,c%q,dt)
  c%outer_seconds(3)=c%outer_seconds(3)+MPI_Wtime()-tick_boundary
  c%history=q0;c%dt_old=dt;c%time=c%time+dt;c%step=c%step+1
  c%phase_seconds(5)=c%phase_seconds(5)+MPI_Wtime()-tick;tick=MPI_Wtime()
  audit=c%phi;call c%faces%sync(audit,st,jump);call demand(st==0,'post CT shared face audit')
  scale=maxval(abs(c%phi));errors(1)=maxval(abs(cube_divergence(c%phi,n)-before))/max(tiny(scale),scale)
  ! sync returns a global jump; normalize the actual local correction by its local flux scale.
  errors(2)=maxval(abs(audit-c%phi))/max(tiny(scale),scale)
  errors(3)=maxval(abs(c%q(1:n(1),1:n(2),1:n(3),:)))
  call demand(all(ieee_is_finite(c%q(1:n(1),1:n(2),1:n(3),:))),'finite updated state')
  call MPI_Allreduce(MPI_IN_PLACE,errors,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,st)
  if(maxval(errors(1:2))>=1.e-11_real64) then
   write(*,*) 'YY_CT_FAILURE rank/step/time/dt/errors/scale/jump=',c%rank,c%step,c%time,dt,errors,scale,jump
   write(failure_name,'("yy_ct_failure_rank",i0.3,".bin")') c%rank
   open(newunit=failure_unit,file=trim(failure_name),status='new',access='stream', &
    form='unformatted',iostat=failure_io)
   if(failure_io==0) then
    write(failure_unit) 'YYCTF001',int([c%rank,c%layout%face,c%step,n,c%layout%first],int32), &
     c%time,dt,errors,scale,jump,c%phi_history,c%phi,emf,audit,before
    close(failure_unit)
   endif
  endif
  call demand(maxval(errors(1:2))<1.e-11_real64,'native CT/shared invariant')
  c%phase_seconds(6)=c%phase_seconds(6)+MPI_Wtime()-tick
  if(clock_reached(c%time,stop))then
   write(*,*)'ACTUAL_PLANE_PACKING_WORD_PASS rank/words=',c%rank,plane_words_compared
   call release_rotation_work()
  endif
  c%phase_seconds(8)=c%phase_seconds(8)+MPI_Wtime()-advance_start
 end subroutine
 subroutine model_dump_mi(c,prefix)
  class(cube_model),intent(in)::c
  character(*),intent(in)::prefix
#ifdef LFM_CUBE_NATIVE_MI
  real(real64),allocatable::data(:,:,:,:)
  character(1024)::name
  integer::u,st
  if(.not.c%mi_enabled.or.c%rank/=0) return
  call c%mi%snapshot(data,st);call demand(st==0,'native MI snapshot state')
  write(name,'(a,"-mi-step",i7.7,".bin")') trim(prefix),c%step
  open(newunit=u,file=trim(name),access='stream',form='unformatted',status='new')
  write(u) 'LFMCUMI1',int(shape(data),int32),int(c%step,int32),c%time,c%mi%last_time,6.5e6_real64,c%mi%sigma,data
  close(u)
#endif
 end subroutine
 subroutine model_profile(c,seconds)
  class(cube_model),intent(in)::c
  real(real64),intent(out)::seconds(11)
  real(real64)::local(11),trap_massmax,trap_energymax
  integer(c_size_t)::trap_calls,trap_faces
  integer::st
  if(c%mode>=1.and.c%layout%coord(1)==0) then
   call yy_flux_trap_statistics(trap_calls,trap_faces,trap_massmax,trap_energymax)
   write(*,*) 'YY_INNER_FLUX_STATS ',c%rank,c%step,c%time,trap_calls,trap_faces,trap_massmax,trap_energymax
  endif
  write(*,*) 'EQUIVALENT_COST_RANK_STEP_TIME_BASE_PRED_FACE_REFERENCE_MAGWALL_BFSWEEP_COROTATION_WALL_ROT_WORK ',c%rank,c%step,c%time,equivalent_cost_s(:8)
  local=[c%phase_seconds,c%io_seconds,c%diagnostic_seconds,c%cleanup_seconds]
  call MPI_Allreduce(local,seconds,11,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,st)
 end subroutine
 subroutine model_diagnose(c,stats)
  class(cube_model),intent(in)::c
  real(real64),intent(out)::stats(6)
  real(real64)::local(6)
  integer::n(3),st
  n=c%n
  local=[minval(c%q(1:n(1),1:n(2),1:n(3),1)),minval(c%q(1:n(1),1:n(2),1:n(3),5)), &
   maxval(abs(c%q(1:n(1),1:n(2),1:n(3),2:4))),maxval(abs(c%q(1:n(1),1:n(2),1:n(3),6:8)+c%background_cell)), &
   maxval(abs(cube_divergence(c%phi,n)))/max(tiny(1._real64),maxval(abs(c%phi))), &
   sum(c%q(1:n(1),1:n(2),1:n(3),1)*c%p%legacy_volume)]
  call MPI_Allreduce(local(1:2),stats(1:2),2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,st)
  call MPI_Allreduce(local(3:5),stats(3:5),3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,st)
  call MPI_Allreduce(local(6),stats(6),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,st)
 end subroutine
 subroutine model_dump(c,prefix)
  class(cube_model),intent(in)::c
  character(*),intent(in)::prefix
  character(1024)::name
  integer::u
  real(real64),allocatable::visible(:,:,:,:)
  write(name,'(a,"-",i0.5,"-r",i3.3,".bin")') trim(prefix),c%step,c%rank
  open(newunit=u,file=trim(name),access='stream',form='unformatted',status='new')
  write(u) 'LFMCUB01',int(c%n,int32),int(c%layout%face,int32),int(c%layout%first,int32),c%time,c%rp
  visible=c%q(1:c%n(1),1:c%n(2),1:c%n(3),:)
  visible(:,:,:,6:8)=visible(:,:,:,6:8)+c%background_cell
  write(u) c%p%x,visible,c%phi
  close(u)
 end subroutine
 ! CPC03 persists only independent evolving state. Geometry, cell B, halos,
 ! CT cuboid padding, held-MI maps/coordinates/conductances and event counters
 ! are reconstructed from immutable inputs by the existing numerical paths.
 integer function checkpoint_group_size(c) result(g)
  class(cube_model),intent(in)::c
  character(32)::value
  integer::ios,lo,hi,st
  g=16;call get_environment_variable('CUBE_CP_GROUP_RANKS',value,status=ios)
  if(ios==0.and.len_trim(value)>0) then
   read(value,*,iostat=ios) g
   call demand(ios==0,'integer checkpoint group policy')
  endif
  call demand(ios/=2.and.g>=1.and.g<=128,'checkpoint group size1..128')
  call MPI_Allreduce(g,lo,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,st)
  call MPI_Allreduce(g,hi,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,st)
  call demand(lo==hi,'checkpoint group policy identical on all ranks')
 end function
 subroutine checkpoint_check(failed,message)
  integer,intent(in)::failed
  character(*),intent(in)::message
  integer::global,st
  call MPI_Allreduce(failed,global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,st)
  call demand(st==MPI_SUCCESS.and.global==0,message)
 end subroutine
 function checkpoint_crc(data) result(crc)
  real(real64),intent(in)::data(:)
  integer(int64)::crc,bits,entry
  integer(int64),save::table(0:255)=0
  logical,save::ready=.false.
  integer::i,j,k,index
  ! CRC64-ECMA, polynomial 0x42F0E1EBA9EA3693, init/xorout0.
  ! Consume each IEEE word least-significant byte first; no integer overflow.
  if(.not.ready) then
   do i=0,255
    entry=shiftl(int(i,int64),56)
    do j=1,8
     if(btest(entry,63)) then
      entry=ieor(shiftl(entry,1),int(z'42F0E1EBA9EA3693',int64))
     else
      entry=shiftl(entry,1)
     endif
    enddo
    table(i)=entry
   enddo
   ready=.true.
  endif
  crc=0
  do i=1,size(data)
   bits=transfer(data(i),bits)
   do k=0,7
    index=int(ieor(shiftr(crc,56),iand(shiftr(bits,8*k),255_int64)))
    crc=ieor(shiftl(crc,8),table(index))
   enddo
  enddo
 end function
 subroutine checkpoint_faces(phi,n,buffer,position,loading)
  integer,intent(in)::n(3)
  real(real64),intent(inout)::phi(0:,0:,0:,:),buffer(:)
  integer,intent(inout)::position
  logical,intent(in)::loading
  integer::d,i,j,k,lo(3),hi(3)
  if(loading) phi=0 ! Unused CT cuboid entries are identically zero.
  do d=1,3
   lo=1;lo(d)=0;hi=n
   do k=lo(3),hi(3)
    do j=lo(2),hi(2)
     do i=lo(1),hi(1)
      position=position+1
      if(loading) then
       phi(i,j,k,d)=buffer(position)
      else
       buffer(position)=phi(i,j,k,d)
      endif
     enddo
    enddo
   enddo
  enddo
 end subroutine
 subroutine checkpoint_pack(c,buffer)
  class(cube_model),intent(inout)::c
  real(real64),allocatable,intent(out)::buffer(:)
  real(real64),allocatable::rebuilt(:,:,:,:),padded(:,:,:,:)
  integer::cells,n(3),p,nf
  n=c%n;cells=product(n)
  nf=(n(1)+1)*n(2)*n(3)+n(1)*(n(2)+1)*n(3)+n(1)*n(2)*(n(3)+1)
  allocate(buffer(10*cells+2*nf),rebuilt(n(1),n(2),n(3),3))
  ! Refuse a format reduction unless the live derived B is exactly reproducible.
  call native_recover(c%p,c%phi-c%background_flux,c%correction,rebuilt)
  call demand(all(transfer(rebuilt,[0_int64],size(rebuilt))== &
   transfer(c%q(1:n(1),1:n(2),1:n(3),6:8),[0_int64],size(rebuilt))),'checkpoint reconstruct current B bytes')
  call native_recover(c%p,c%phi_history-c%background_flux,c%correction,rebuilt)
  call demand(all(transfer(rebuilt,[0_int64],size(rebuilt))== &
   transfer(c%history(1:n(1),1:n(2),1:n(3),6:8),[0_int64],size(rebuilt))),'checkpoint reconstruct history B bytes')
  ! Verify even unused zero padding rather than silently discarding information.
  padded=c%phi;p=10*cells;buffer=0
  buffer(1:5*cells)=reshape(c%q(1:n(1),1:n(2),1:n(3),1:5),[5*cells])
  buffer(5*cells+1:10*cells)=reshape(c%history(1:n(1),1:n(2),1:n(3),1:5),[5*cells])
  call checkpoint_faces(c%phi,n,buffer,p,.false.)
  call checkpoint_faces(c%phi_history,n,buffer,p,.false.)
  call demand(p==size(buffer),'compact checkpoint field accounting')
  p=10*cells;call checkpoint_faces(padded,n,buffer,p,.true.)
  call demand(all(transfer(padded,[0_int64],size(padded))== &
   transfer(c%phi,[0_int64],size(padded))),'checkpoint CT padding reconstructs exactly')
  call checkpoint_faces(padded,n,buffer,p,.true.)
  call demand(all(transfer(padded,[0_int64],size(padded))== &
   transfer(c%phi_history,[0_int64],size(padded))),'checkpoint old CT padding reconstructs exactly')
 end subroutine
 subroutine checkpoint_unpack(c,buffer)
  class(cube_model),intent(inout)::c
  real(real64),intent(inout)::buffer(:)
  integer::n(3),cells,p
  n=c%n;cells=product(n);c%q=0;c%history=0
  ! Valid temporary ghost placeholders; both full halos are rebuilt below.
  c%q(:,:,:,1)=c%sw(1);c%q(:,:,:,5)=c%sw(5)
  c%history(:,:,:,1)=c%sw(1);c%history(:,:,:,5)=c%sw(5)
  c%q(1:n(1),1:n(2),1:n(3),1:5)=reshape(buffer(:5*cells),[n,5])
  c%history(1:n(1),1:n(2),1:n(3),1:5)=reshape(buffer(5*cells+1:10*cells),[n,5])
  p=10*cells;call checkpoint_faces(c%phi,n,buffer,p,.true.)
  call checkpoint_faces(c%phi_history,n,buffer,p,.true.)
  call demand(p==size(buffer),'restored compact field accounting')
  call native_recover(c%p,c%phi-c%background_flux,c%correction,c%q(1:n(1),1:n(2),1:n(3),6:8))
  call native_recover(c%p,c%phi_history-c%background_flux,c%correction,c%history(1:n(1),1:n(2),1:n(3),6:8))
  ! Only active history participates in the predictor; predictor halos are
  ! overwritten by existing guards before HYDRO/BFIELD read them.
  call c%guards(c%history);call c%guards(c%q)
 end subroutine
 subroutine checkpoint_mi_clock(c)
  class(cube_model),intent(inout)::c
  c%mi_next=0;c%mi_updates=0
  do while(clock_reached(c%time,c%mi_next))
   c%mi_next=c%mi_next+c%mi_cadence;c%mi_updates=c%mi_updates+1
  enddo
 end subroutine
 subroutine model_checkpoint(c,prefix,fingerprint)
  class(cube_model),intent(inout)::c
  character(*),intent(in)::prefix,fingerprint
  type(MPI_Comm)::group
  real(real64),allocatable::buffer(:),received(:),mi_state(:)
  real(real64)::clock(2),shared(2),next_mi
  integer(int64),allocatable::crcs(:)
  integer(int64)::crc,mi_crc,step,words,file_bytes
  character(1024)::name,stem
  integer::u,io,st,failed,g,groups,gid,grank,members,first,i,updates
  logical::exists
  g=checkpoint_group_size(c);groups=(c%ranks+g-1)/g;gid=c%rank/g;first=gid*g
  call checkpoint_pack(c,buffer);words=int(size(buffer),int64)
  call demand(words*g<huge(1),'checkpoint bounded MPI count')
  clock=[c%time,c%dt_old];shared=clock
  call MPI_Bcast(shared,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,st)
  call demand(all(transfer(clock,[0_int64],2)==transfer(shared,[0_int64],2)),'checkpoint shared exact clock')
  step=int(c%step,int64)
  write(stem,'(a,"-step",i0.7)') trim(prefix),step
  write(name,'(a,"-g",i0.6,".cp")') trim(stem),gid
  failed=0;inquire(file=trim(stem)//'.complete',exist=exists);if(exists) failed=1
  call checkpoint_check(failed,'checkpoint refuses existing complete epoch')
  call MPI_Comm_split(MPI_COMM_WORLD,gid,c%rank,group,st)
  call MPI_Comm_rank(group,grank,st);call MPI_Comm_size(group,members,st)
  allocate(received(merge(size(buffer)*members,1,grank==0)))
  call MPI_Gather(buffer,size(buffer),MPI_DOUBLE_PRECISION,received,size(buffer),MPI_DOUBLE_PRECISION,0,group,st)
  failed=merge(0,1,st==MPI_SUCCESS);crc=0
  if(grank==0.and.failed==0) then
   crc=checkpoint_crc([clock,received])
   open(newunit=u,file=trim(name)//'.tmp',access='stream',form='unformatted',status='new',iostat=io)
   if(io==0) then
    write(u,iostat=io) 'LFMCPS03',fingerprint,step,clock, &
     int([c%mode,c%ranks,g,gid,first,members,c%n,c%layout%global,c%layout%dims],int32),words,received,crc
    if(io/=0) failed=1
    close(u,iostat=io)
   endif
   if(io/=0) failed=1
  endif
  call checkpoint_check(failed,'checkpoint shard temporary write')
  call MPI_Bcast(crc,1,MPI_INTEGER8,0,group,st)
  allocate(crcs(c%ranks));call MPI_Gather(crc,1,MPI_INTEGER8,crcs,1,MPI_INTEGER8,0,MPI_COMM_WORLD,st)
  mi_crc=0;allocate(mi_state(0))
#ifdef LFM_CUBE_NATIVE_MI
  if(c%mi_enabled) then
   deallocate(mi_state);call c%mi%checkpoint_export(mi_state,st)
   call demand(st==0,'minimal held MI export')
   next_mi=c%mi_next;updates=c%mi_updates
   call checkpoint_mi_clock(c)
   call demand(c%mi_next==next_mi.and.c%mi_updates==updates,'MI schedule is exactly reconstructable')
   mi_crc=checkpoint_crc(mi_state)
   call MPI_Bcast(mi_crc,1,MPI_INTEGER8,0,MPI_COMM_WORLD,st)
   call demand(checkpoint_crc(mi_state)==mi_crc,'shared independent MI state')
   if(c%rank==0) then
    open(newunit=u,file=trim(stem)//'.mi.tmp',access='stream',form='unformatted',status='new',iostat=io)
    if(io==0) then
     write(u,iostat=io) 'LFMCPI03',fingerprint,step,c%time,int(size(mi_state),int64),mi_state,mi_crc
     if(io/=0) failed=1
     close(u,iostat=io)
    endif
    if(io/=0) failed=1
   endif
  endif
#endif
  call checkpoint_check(failed,'checkpoint single-owner MI temporary write')
  if(grank==0) then
   inquire(file=trim(name),exist=exists)
   if(exists) then
    failed=1
   else
    call rename(trim(name)//'.tmp',trim(name),io);if(io/=0) failed=1
   endif
  endif
  if(c%rank==0.and.c%mi_enabled) then
   inquire(file=trim(stem)//'.mi',exist=exists)
   if(exists) then
    failed=1
   else
    call rename(trim(stem)//'.mi.tmp',trim(stem)//'.mi',io);if(io/=0) failed=1
   endif
  endif
  call checkpoint_check(failed,'checkpoint shards committed before manifest')
  if(c%rank==0) then
   open(newunit=u,file=trim(stem)//'.complete.tmp',status='new',iostat=io)
   if(io==0) then
    write(u,'(a)',iostat=io) fingerprint
    if(io==0) write(u,*,iostat=io) step,c%time,c%ranks
    if(io==0) write(u,*,iostat=io) 'LFMCPC03',g,groups,c%mode
    if(io==0) write(u,*,iostat=io) c%layout%global,c%layout%dims
    if(io==0) write(u,*,iostat=io) size(mi_state),mi_crc
    do i=0,groups-1
     file_bytes=172_int64+8*words*min(g,c%ranks-i*g)
     if(io==0) write(u,*,iostat=io) i,i*g,min(g,c%ranks-i*g),file_bytes,crcs(i*g+1)
    enddo
    if(io/=0) failed=1
    close(u,iostat=io)
   endif
   if(io/=0) failed=1
   if(failed==0) call rename(trim(stem)//'.complete.tmp',trim(stem)//'.complete',io)
   if(io/=0) failed=1
  endif
  call checkpoint_check(failed,'checkpoint final complete manifest publication')
  call MPI_Comm_free(group,st)
 end subroutine
 subroutine model_restore(c,stem,fingerprint)
  class(cube_model),intent(inout)::c
  character(*),intent(in)::stem,fingerprint
  type(MPI_Comm)::group
  character(1024)::name
  character(64)::stored
  character(8)::magic
  real(real64),allocatable::buffer(:),received(:),mi_state(:)
  real(real64)::epoch,clock(2)
  integer(int32)::meta(15)
  integer(int64),allocatable::checks(:),bytes(:)
  integer(int64)::step,diskstep,words,actual_words,crc,mi_crc,size_bytes
  integer::u,io,st,failed,g,groups,gid,grank,members,first,mode,ranks,global(3),dims(3),mi_count,i,j,k
  integer::fields(3),n(3),cells,p,extra
  failed=0;g=0
  ! Root reads the small manifest once, then broadcasts metadata.
  if(c%rank==0) then
   open(newunit=u,file=trim(stem)//'.complete',status='old',iostat=io)
   if(io==0) then
    read(u,'(a)',iostat=io) stored
    if(io==0) read(u,*,iostat=io) step,epoch,ranks
    if(io==0) read(u,*,iostat=io) magic,g,groups,mode
    if(io<0) g=-1 ! Explicit old-format reader for controlled migration/proofs.
    if(io>0) failed=1
    close(u)
   else
    failed=1
   endif
  endif
  call checkpoint_check(failed,'checkpoint readable complete manifest required')
  call MPI_Bcast(g,1,MPI_INTEGER,0,MPI_COMM_WORLD,st)
  if(g==-1) then
   call model_restore_legacy(c,stem,fingerprint);return
  endif
  call demand(g>=1.and.g<=128,'checkpoint group policy range')
  groups=(c%ranks+g-1)/g;allocate(checks(groups),bytes(groups))
  if(c%rank==0) then
   open(newunit=u,file=trim(stem)//'.complete',status='old',iostat=io)
   if(io==0) read(u,'(a)',iostat=io) stored
   if(io==0) read(u,*,iostat=io) step,epoch,ranks
   if(io==0) read(u,*,iostat=io) magic,g,i,mode
   if(io==0) then
    if(stored/=fingerprint.or.magic/='LFMCPC03'.or.i/=groups.or.ranks/=c%ranks.or.mode/=c%mode) failed=1
   endif
   if(io==0) read(u,*,iostat=io) global,dims
   if(io==0) then
    if(any(global/=c%layout%global).or.any(dims/=c%layout%dims)) failed=1
   endif
   if(io==0) read(u,*,iostat=io) mi_count,mi_crc
   do i=1,groups
    if(io/=0) exit
    read(u,*,iostat=io) fields,bytes(i),checks(i)
    if(io==0) then
     if(any(fields/=[i-1,(i-1)*g,min(g,c%ranks-(i-1)*g)])) failed=1
    endif
   enddo
   if(io/=0) failed=1
   if(io==0) then
    read(u,*,iostat=io) extra
    if(io>=0) failed=1
   endif
   close(u)
  endif
  call checkpoint_check(failed,'checkpoint manifest schema/layout/ownership')
  call MPI_Bcast(step,1,MPI_INTEGER8,0,MPI_COMM_WORLD,st)
  call MPI_Bcast(epoch,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,st)
  call MPI_Bcast(mi_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,st)
  call MPI_Bcast(mi_crc,1,MPI_INTEGER8,0,MPI_COMM_WORLD,st)
  call MPI_Bcast(checks,groups,MPI_INTEGER8,0,MPI_COMM_WORLD,st)
  call MPI_Bcast(bytes,groups,MPI_INTEGER8,0,MPI_COMM_WORLD,st)
  call demand(step>0.and.step<=huge(c%step).and.epoch>0,'checkpoint supported exact step/epoch')
  n=c%n;cells=product(n)
  words=10_int64*cells+2_int64*((n(1)+1)*n(2)*n(3)+n(1)*(n(2)+1)*n(3)+n(1)*n(2)*(n(3)+1))
  gid=c%rank/g;first=gid*g
  call MPI_Comm_split(MPI_COMM_WORLD,gid,c%rank,group,st)
  call MPI_Comm_rank(group,grank,st);call MPI_Comm_size(group,members,st)
  allocate(buffer(words),received(merge(int(words)*members,1,grank==0)))
  write(name,'(a,"-g",i0.6,".cp")') trim(stem),gid
  clock=0
  if(grank==0) then
   inquire(file=trim(name),size=size_bytes,iostat=io)
   if(io/=0.or.size_bytes/=bytes(gid+1).or.size_bytes/=172+8*words*members) failed=1
   if(failed==0) then
    open(newunit=u,file=trim(name),access='stream',form='unformatted',status='old',iostat=io)
    if(io==0) then
     read(u,iostat=io) magic,stored,diskstep,clock,meta,actual_words,received,crc
     if(io==0) then
      if(magic/='LFMCPS03'.or.stored/=fingerprint.or.diskstep/=step.or.clock(1)/=epoch.or.clock(2)<=0) failed=1
      if(any(meta/=int([c%mode,c%ranks,g,gid,first,members,n,c%layout%global,c%layout%dims],int32))) failed=1
      if(actual_words/=words.or.crc/=checks(gid+1).or.checkpoint_crc([clock,received])/=crc) failed=1
     endif
     close(u)
    endif
    if(io/=0) failed=1
   endif
  endif
  call checkpoint_check(failed,'checkpoint shard length/schema/CRC64 before state install')
  call MPI_Bcast(clock,2,MPI_DOUBLE_PRECISION,0,group,st)
  call MPI_Scatter(received,int(words),MPI_DOUBLE_PRECISION,buffer,int(words),MPI_DOUBLE_PRECISION,0,group,st)
  if(st/=MPI_SUCCESS.or..not.all(ieee_is_finite(buffer))) failed=1
  if(minval(buffer(:cells))<=0.or.minval(buffer(4*cells+1:5*cells))<=0) failed=1
  if(minval(buffer(5*cells+1:6*cells))<=0.or.minval(buffer(9*cells+1:10*cells))<=0) failed=1
  call checkpoint_check(failed,'checkpoint finite positive independent fluid/CT state')
  call demand(mi_count>=0.and.mi_count<1000000,'checkpoint independent MI count bound')
  allocate(mi_state(mi_count))
  if(c%rank==0.and.c%mi_enabled) then
   inquire(file=trim(stem)//'.mi',size=size_bytes,iostat=io)
   if(io/=0.or.size_bytes/=104_int64+8_int64*mi_count) failed=1
   if(failed==0) then
    open(newunit=u,file=trim(stem)//'.mi',access='stream',form='unformatted',status='old',iostat=io)
    if(io==0) then
     read(u,iostat=io) magic,stored,diskstep,epoch,actual_words,mi_state,crc
     if(io==0) then
      if(magic/='LFMCPI03'.or.stored/=fingerprint.or.diskstep/=step.or.epoch/=clock(1)) failed=1
      if(actual_words/=mi_count.or.crc/=mi_crc.or.checkpoint_crc(mi_state)/=crc) failed=1
     endif
     close(u)
    endif
    if(io/=0) failed=1
   endif
  endif
  if(.not.c%mi_enabled.and.mi_count/=0) failed=1
  call checkpoint_check(failed,'checkpoint owner MI length/schema/CRC64 before install')
  if(mi_count>0) call MPI_Bcast(mi_state,mi_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,st)
  c%step=int(step);c%time=clock(1);c%dt_old=clock(2)
#ifdef LFM_CUBE_NATIVE_MI
  if(c%mi_enabled) then
   call c%mi%checkpoint_import(mi_state,st);call demand(st==0,'rebuild native held MI from independent state')
   call demand(c%mi%last_step<=c%step.and.c%mi%last_time<=c%time,'restored held MI epoch')
   call checkpoint_mi_clock(c)
  endif
#endif
  call checkpoint_unpack(c,buffer)
  call MPI_Comm_free(group,st)
 end subroutine

 subroutine model_restore_legacy(c,stem,fingerprint)
  class(cube_model),intent(inout)::c
  character(*),intent(in)::stem,fingerprint
  character(1024)::name
  character(8)::magic
  character(64)::stored,marker
  integer(int32)::n(3),global(3),first(3),face,ranks,step,mode,extra
  real(real64)::time,dt,rp,moment,sw(8),low,high
  real(real64)::shared_clock(3)
  integer::u,st
  logical::complete
#ifdef LFM_CUBE_NATIVE_MI
  real(real64),allocatable::mi_state(:)
  integer(int32)::mi_count,mi_updates,shared_mi_count
  real(real64)::mi_next,mi_cadence
  real(real64)::mi_clock(3)
  integer(int64),allocatable::mi_bytes(:),shared_mi_bytes(:)
#endif
  inquire(file=trim(stem)//'.complete',exist=complete);call demand(complete,'complete checkpoint marker required')
  open(newunit=u,file=trim(stem)//'.complete',status='old');read(u,'(a)') marker;close(u)
  call demand(marker==fingerprint,'checkpoint marker provenance mismatch')
  write(name,'(a,"-r",i3.3,".cp")') trim(stem),c%rank
  open(newunit=u,file=trim(name),access='stream',form='unformatted',status='old')
  read(u,iostat=st) magic,stored,n,global,first,face,ranks,step,mode,time,dt,rp,moment,sw
  call demand(st==0,'checkpoint header complete')
  call demand(magic==merge('LFMCPC02','LFMCPC01',c%mi_enabled).and.stored==fingerprint,'checkpoint version/provenance')
  call demand(all(n==c%n).and.all(global==c%layout%global).and.all(first==c%layout%first) &
    .and.face==c%layout%face.and.ranks==c%ranks,'checkpoint layout mismatch')
  call demand(mode==c%mode.and.rp==c%rp.and.moment==c%moment.and.all(sw==c%sw),'checkpoint physics mismatch')
  read(u,iostat=st) c%q,c%history,c%phi,c%phi_history
  call demand(st==0,'checkpoint state complete')
#ifdef LFM_CUBE_NATIVE_MI
  if(c%mi_enabled) then
   read(u,iostat=st) mi_count,mi_next,mi_cadence,mi_updates
   call demand(st==0,'native MI checkpoint header')
   call demand(mi_count>0.and.mi_count<1000000.and.mi_cadence==10._real64.and.mi_next>time &
    .and.mi_updates>0,'native MI checkpoint clocks/shape')
   shared_mi_count=mi_count
   call MPI_Bcast(shared_mi_count,1,MPI_INTEGER4,0,MPI_COMM_WORLD,st)
   call demand(shared_mi_count==mi_count,'shared MI checkpoint payload size')
   allocate(mi_state(mi_count));read(u,iostat=st) mi_state;call demand(st==0,'native MI checkpoint payload')
   ! All ranks must hold the same scheduling and solver state. Otherwise
   ! one rank can skip a collective that another enters after restart.
   mi_clock=[mi_next,mi_cadence,real(mi_updates,real64)]
   call MPI_Bcast(mi_clock,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,st)
   call demand(all(mi_clock==[mi_next,mi_cadence,real(mi_updates,real64)]),'shared MI checkpoint scheduler')
   mi_bytes=transfer(mi_state,[0_int64],size(mi_state));shared_mi_bytes=mi_bytes
   call MPI_Bcast(shared_mi_bytes,size(shared_mi_bytes),MPI_INTEGER8,0,MPI_COMM_WORLD,st)
   call demand(all(mi_bytes==shared_mi_bytes),'shared MI checkpoint state bytes')
   call c%mi%import(mi_state,st);call demand(st==0,'native MI warm-start/held-state restore')
   call demand(c%mi%last_time<=time.and.c%mi%last_step<=step,'native MI checkpoint solution epoch')
   c%mi_next=mi_next;c%mi_cadence=mi_cadence;c%mi_updates=mi_updates
  endif
#endif
  read(u,iostat=st) extra;call demand(st<0,'unexpected checkpoint payload');close(u)
  call demand(all(ieee_is_finite(c%q)).and.all(ieee_is_finite(c%history)) &
   .and.all(ieee_is_finite(c%phi)).and.all(ieee_is_finite(c%phi_history)),'checkpoint finite state')
  call demand(minval(c%q(:,:,:,1))>0.and.minval(c%q(:,:,:,5))>0 &
    .and.minval(c%history(:,:,:,1))>0.and.minval(c%history(:,:,:,5))>0,'checkpoint positive state')
  call MPI_Allreduce(time,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,st)
  call MPI_Allreduce(time,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,st)
  call demand(low==high.and.dt>0,'checkpoint shared clock')
  shared_clock=[time,dt,real(step,real64)]
  call MPI_Bcast(shared_clock,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,st)
  call demand(all(shared_clock==[time,dt,real(step,real64)]),'checkpoint shared step/history clock')
  c%step=step;c%time=time;c%dt_old=dt
 end subroutine
 ! Outer boundary adapted from kaiju-master wind.F90 (WindBC, WindEFix,
 ! NudgeSW). Keep the native wind clock and the CT-owned physical face flux.
 ! Ghost recovery is the existing calibrated-BZZ algebra: inverse of the
 ! three averaged area vectors, as in cube_bzz_calibration/native_recover.
 function kaiju_wind_weight(sw,x) result(w)
  real(real64),intent(in)::sw(8),x(3)
  real(real64)::w,angle,s,z
  z=-dot_product(sw(2:4),x)/(norm2(sw(2:4))*norm2(x))
  angle=acos(max(-1._real64,min(1._real64,z)))*180._real64/acos(-1._real64)
  w=1
  if(angle>=180._real64) then
   w=0
  else if(angle>90._real64) then
   s=(angle-90._real64)/90._real64
   w=exp(1._real64-1._real64/(1._real64-s*s))
  endif
 end function
 function kaiju_outer_state(sw,time,x,donor,normal,weight,switch_s,bz_after,reference_offset_s) result(q)
  real(real64),intent(in)::sw(8),time,x(3),donor(8),normal(3)
  real(real64),optional,intent(in)::weight
  real(real64),optional,intent(in)::switch_s,bz_after,reference_offset_s
  real(real64)::q(8),wind(8),outflow(8),w,fast,va,vn,pressure
  wind=earth_wind(sw,time,x,switch_s,bz_after,reference_offset_s)
  if(present(weight)) then
   w=weight
  else
   w=kaiju_wind_weight(sw,x)
  endif
  outflow=donor
  va=norm2(donor(6:8))/sqrt(4*acos(-1._real64)*donor(1))
  va=ca*va/sqrt(ca*ca+va*va);fast=sqrt(donor(5)**2+va*va)
  vn=dot_product(donor(2:4),normal)
  if(vn<fast) outflow(2:4)=outflow(2:4)+(fast-vn)*normal
  q=w*wind+(1-w)*outflow
  pressure=w*wind(1)*wind(5)**2/gamma+(1-w)*donor(1)*donor(5)**2/gamma
  call demand(q(1)>0.and.pressure>0,'positive Kaiju mixed boundary state')
  q(5)=sqrt(gamma*pressure/q(1))
 end function
 subroutine prepare_outer_boundary(c)
  class(cube_model),intent(inout)::c
  integer::n(3),r,i,j,k,d,h,idx(3),hi(3),low(3),rmax,jmax,kmax,st,imin,jmin,kmin,perp(3)
  real(real64)::v(3,4),matrix(3,3),inverse(3,3),a0(3),b0(3),edge(3),x(3),length1,length2
  if(c%mode<2.or.c%layout%coord(1)/=c%layout%dims(1)-1) return
  if(allocated(c%outer_area)) return
  n=c%n
  allocate(c%outer_area(3,0:ng+1,1-ng:n(2)+ng+1,1-ng:n(3)+ng+1,3), &
   c%outer_xyz(3,0:ng+1,1-ng:n(2)+ng+1,1-ng:n(3)+ng+1,3), &
   c%outer_inverse(3,3,ng,1-ng:n(2)+ng,1-ng:n(3)+ng), &
   c%outer_curl_inverse(3,3,n(1):n(1)+1,0:n(2)+1,0:n(3)+1), &
   c%outer_weight(0:ng+1,1-ng:n(2)+ng+1,1-ng:n(3)+ng+1,3), &
   c%outer_magnitude(0:ng+1,1-ng:n(2)+ng+1,1-ng:n(3)+ng+1,3), &
   c%outer_cell_xyz(3,0:ng,1-ng:n(2)+ng,1-ng:n(3)+ng), &
   c%outer_cell_weight(0:ng,1-ng:n(2)+ng,1-ng:n(3)+ng), &
   c%outer_normal(3,1-ng:n(2)+ng,1-ng:n(3)+ng), &
   c%outer_edge_xyz(3,n(1)-3:n(1),0:n(2),0:n(3),3), &
   c%outer_edge_vector(3,n(1)-3:n(1),0:n(2),0:n(3),3), &
   c%outer_edge_weight(n(1)-3:n(1),0:n(2),0:n(3),3), &
   c%outer_edge_dl(n(1)-3:n(1),0:n(2),0:n(3),3))
  c%outer_area=0;c%outer_xyz=0
  c%outer_weight=0;c%outer_magnitude=0;c%outer_edge_dl=0
  do d=1,3
   rmax=ng;jmax=n(2)+ng;kmax=n(3)+ng
   if(d==1) rmax=ng+1
   if(d==2) jmax=jmax+1
   if(d==3) kmax=kmax+1
   do k=1-ng,kmax
    do j=1-ng,jmax
     do r=0,rmax
      idx=[n(1)+r,j,k]
      if(r==0) idx(1)=n(1)+merge(1,0,d==1)
      call get_face(c,d,idx,v)
      c%outer_area(:,r,j,k,d)=cube_quad(v)
      c%outer_xyz(:,r,j,k,d)=sum(v,dim=2)/4
      c%outer_weight(r,j,k,d)=kaiju_wind_weight(c%sw,c%outer_xyz(:,r,j,k,d))
      c%outer_magnitude(r,j,k,d)=norm2(c%outer_area(:,r,j,k,d))
      call demand(norm2(c%outer_area(:,r,j,k,d))>0,'positive boundary face area')
      call demand(norm2(c%outer_xyz(:,r,j,k,d))>8*c%rp,'boundary outside compact B0')
     enddo
    enddo
   enddo
  enddo
  do k=1-ng,n(3)+ng
   do j=1-ng,n(2)+ng
    do r=1,ng
     do d=1,3
      hi=[r,j,k];hi(d)=hi(d)+1
      matrix(d,:)=(c%outer_area(:,r,j,k,d)+c%outer_area(:,hi(1),hi(2),hi(3),d))/2
     enddo
     call cube_inverse3(matrix,inverse,st)
     call demand(st==0,'native calibrated BZZ ghost geometry inverse')
     c%outer_inverse(:,:,r,j,k)=inverse
    enddo
   enddo
  enddo
  do k=0,n(3)+1
   do j=0,n(2)+1
    do i=n(1),n(1)+1
     do d=1,3
      low=[i,j,k];hi=low;low(d)=low(d)-1;hi(d)=hi(d)+1
      matrix(d,:)=boundary_cell_center(c,hi)-boundary_cell_center(c,low)
     enddo
     call cube_inverse3(matrix,inverse,st)
     call demand(st==0,'boundary Cartesian derivative geometry inverse')
     c%outer_curl_inverse(:,:,i,j,k)=inverse
    enddo
   enddo
  enddo
  do k=1-ng,n(3)+ng
   do j=1-ng,n(2)+ng
    c%outer_normal(:,j,k)=c%area(:,n(1)+1,j,k,1)/c%area_magnitude(n(1)+1,j,k,1)
    do r=0,ng
     i=n(1)+r
     x=(c%xv(:,i-1,j-1,k-1)+c%xv(:,i,j-1,k-1)+c%xv(:,i-1,j,k-1)+c%xv(:,i,j,k-1) &
       +c%xv(:,i-1,j-1,k)+c%xv(:,i,j-1,k)+c%xv(:,i-1,j,k)+c%xv(:,i,j,k))/8
     c%outer_cell_xyz(:,r,j,k)=x
     c%outer_cell_weight(r,j,k)=kaiju_wind_weight(c%sw,x)
    enddo
   enddo
  enddo
  c%outer_wind_speed=norm2(c%sw(2:4))
  do d=1,3
   imin=n(1)-3;if(d==1) imin=n(1)-2
   jmin=0;if(d==2) jmin=1
   kmin=0;if(d==3) kmin=1
   do k=kmin,n(3)
    do j=jmin,n(2)
     do i=imin,n(1)
      idx=[i,j,k];low=idx;low(d)=low(d)-1
      a0=c%xv(:,low(1),low(2),low(3));b0=c%xv(:,i,j,k);edge=b0-a0;x=(a0+b0)/2
      c%outer_edge_xyz(:,i,j,k,d)=x;c%outer_edge_vector(:,i,j,k,d)=edge
      c%outer_edge_weight(i,j,k,d)=kaiju_wind_weight(c%sw,x)
      if(i/=n(1).or.d==1) cycle
      perp=low;perp(5-d)=perp(5-d)+1
      length1=norm2(edge);length2=norm2(c%xv(:,perp(1),perp(2),perp(3))-a0)
      c%outer_edge_dl(i,j,k,d)=sqrt(length1*length2)
     enddo
    enddo
   enddo
  enddo
 end subroutine
 function boundary_cell_center(c,p) result(x)
  class(cube_model),intent(in)::c
  integer,intent(in)::p(3)
  real(real64)::x(3)
  integer::h
  x=0
  do h=0,7
   x=x+c%xv(:,p(1)-1+mod(h,2),p(2)-1+mod(h/2,2),p(3)-1+h/4)
  enddo
  x=x/8
 end function
 function outer_donor_flux(c,q,phi,d,j,k) result(value)
  class(cube_model),intent(in)::c
  real(real64),intent(in)::q(1-ng:,1-ng:,1-ng:,:),phi(0:,0:,0:,:)
  integer,intent(in)::d,j,k
  real(real64)::value,b(3)
  integer::n(3),a(3),z(3)
  n=c%n
  select case(d)
   case(1)
    if(j>=1.and.j<=n(2).and.k>=1.and.k<=n(3)) then
     value=phi(n(1),j,k,1);return
    endif
   case(2)
    if(j>=1.and.j<=n(2)+1.and.k>=1.and.k<=n(3)) then
     value=phi(n(1),j-1,k,2);return
    endif
   case(3)
    if(j>=1.and.j<=n(2).and.k>=1.and.k<=n(3)+1) then
     value=phi(n(1),j,k-1,3);return
    endif
  end select
  ! Angular guards use the already qualified Cartesian halo interpolation.
  a=[n(1),j,k];z=a
  if(d/=1) a(d)=a(d)-1
  a=max([1-ng,1-ng,1-ng],min(n+ng,a));z=max([1-ng,1-ng,1-ng],min(n+ng,z))
  b=(q(a(1),a(2),a(3),6:8)+q(z(1),z(2),z(3),6:8))/2
  value=dot_product(b,c%outer_area(:,0,j,k,d))
 end function
 subroutine outer_face_fields(c,q,phi,bf)
  class(cube_model),intent(inout)::c
  real(real64),intent(inout)::q(1-ng:,1-ng:,1-ng:,:),bf(1-ng:,1-ng:,1-ng:,:)
  real(real64),intent(in)::phi(0:,0:,0:,:)
  real(real64),allocatable::flux(:,:,:,:)
  real(real64)::wind(8),w,donor,source(3),b(3),area
  integer::n(3),r,j,k,d,hi(3),rmax,jmax,kmax
  if(c%mode<2.or.c%layout%coord(1)/=c%layout%dims(1)-1) return
  n=c%n;call demand(allocated(c%outer_area),'boundary geometry cache ready')
  allocate(flux(ng+1,1-ng:n(2)+ng+1,1-ng:n(3)+ng+1,3));flux=0
  do d=1,3
   rmax=ng;jmax=n(2)+ng;kmax=n(3)+ng
   if(d==1) rmax=rmax+1
   if(d==2) jmax=jmax+1
   if(d==3) kmax=kmax+1
   do k=1-ng,kmax
    do j=1-ng,jmax
     donor=outer_donor_flux(c,q,phi,d,j,k)
     do r=1,rmax
      w=c%outer_weight(r,j,k,d)
      wind=earth_wind(c%sw,c%time,c%outer_xyz(:,r,j,k,d), &
       c%wind_switch_s,c%wind_bz_after,c%wind_reference_offset_s)
      flux(r,j,k,d)=w*dot_product(wind(6:8),c%outer_area(:,r,j,k,d))+(1-w)*donor
      if(d==1.and.r==1) flux(r,j,k,d)=donor ! Physical Fn stays CT-owned.
      if(r<=ng.and.j<=n(2)+ng.and.k<=n(3)+ng) then
       area=c%outer_magnitude(r,j,k,d)
       bf(n(1)+r,j,k,d)=flux(r,j,k,d)/area
      endif
     enddo
    enddo
   enddo
  enddo
  do k=1-ng,n(3)+ng
   do j=1-ng,n(2)+ng
    do r=1,ng
     do d=1,3
      hi=[r,j,k];hi(d)=hi(d)+1
      source(d)=(flux(r,j,k,d)+flux(hi(1),hi(2),hi(3),d))/2
     enddo
     b=matmul(c%outer_inverse(:,:,r,j,k),source)
     q(n(1)+r,j,k,6:8)=b
    enddo
   enddo
  enddo
  deallocate(flux)
 end subroutine
 subroutine boundary_cartesian_curl(c,q,curl)
  ! Metric adaptation of Kaiju DiffuseOuter current. Cartesian B differences
  ! along three actual opposite-neighbour centre chords are transformed by
  ! the cached native cube_inverse3. This removes rotating-face-normal terms
  ! and is exact for affine Cartesian B, including constant B, on this mesh.
  class(cube_model),intent(in)::c
  real(real64),intent(in)::q(1-ng:,1-ng:,1-ng:,:)
  real(real64),allocatable,intent(out)::curl(:,:,:,:)
  integer::n(3),i,j,k,d,p(3),lo(3),hi(3)
  real(real64)::difference(3,3),gradient(3,3)
  n=c%n;allocate(curl(3,n(1):n(1)+1,0:n(2)+1,0:n(3)+1));curl=0
  do k=0,n(3)+1
   do j=0,n(2)+1
    do i=n(1),n(1)+1
     p=[i,j,k]
     do d=1,3
      lo=p;lo(d)=lo(d)-1;hi=p;hi(d)=hi(d)+1
      difference(d,:)=q(hi(1),hi(2),hi(3),6:8)-q(lo(1),lo(2),lo(3),6:8)
     enddo
     gradient=matmul(c%outer_curl_inverse(:,:,i,j,k),difference)
     curl(:,i,j,k)=[gradient(2,3)-gradient(3,2),gradient(3,1)-gradient(1,3),gradient(1,2)-gradient(2,1)]
    enddo
   enddo
  enddo
 end subroutine
 subroutine kaiju_outer_emf(c,q,dt,emf,diffusion_only)
  class(cube_model),intent(in)::c
  real(real64),intent(in)::q(1-ng:,1-ng:,1-ng:,:),dt
  real(real64),intent(inout)::emf(0:,0:,0:,:)
  logical,optional,intent(in)::diffusion_only
  real(real64),allocatable::curl(:,:,:,:)
  real(real64)::a0(3),b0(3),edge(3),x(3),wind(8),w,taper,dl,vd,jedge(3),length1,length2
  integer::n(3),d,i,j,k,p(3),lo(3),perp(3),level,gI,global_nr,imin,jmin,kmin
  logical::only
  if(c%mode<2.or.c%layout%coord(1)/=c%layout%dims(1)-1) return
  only=.false.;if(present(diffusion_only)) only=diffusion_only
  n=c%n;global_nr=c%layout%global(1)
  call boundary_cartesian_curl(c,q,curl)
  do d=1,3
   imin=n(1)-3;if(d==1) imin=n(1)-2
   jmin=0;if(d==2) jmin=1
   kmin=0;if(d==3) kmin=1
   do k=kmin,n(3)
    do j=jmin,n(2)
     do i=imin,n(1)
      p=[i,j,k];lo=p;lo(d)=lo(d)-1
      x=c%outer_edge_xyz(:,i,j,k,d);edge=c%outer_edge_vector(:,i,j,k,d)
      w=c%outer_edge_weight(i,j,k,d);gI=c%layout%first(1)+i
      level=global_nr-gI;if(d==1) level=level+1
      call demand(level>=0.and.level<=3,'Kaiju four radial edge levels')
      if(.not.only) then
       taper=w/real(level+1,real64)
       wind=earth_wind(c%sw,c%time,x,c%wind_switch_s,c%wind_bz_after,c%wind_reference_offset_s)
       emf(i,j,k,d)=(1-taper)*emf(i,j,k,d)+taper*dot_product(cube_cross(wind(2:4),wind(6:8)),edge)
      endif
      if(i/=n(1).or.d==1) cycle
      ! Kaiju eta=Vd*dl, dl=sqrt(two tangent edge lengths at one vertex).
      dl=c%outer_edge_dl(i,j,k,d);vd=min(c%outer_wind_speed,ca,.3_real64*dl/dt)
      if(d==2) then
       jedge=(curl(:,n(1),j,k)+curl(:,n(1)+1,j,k)+curl(:,n(1),j,k+1)+curl(:,n(1)+1,j,k+1))/4
      else
       jedge=(curl(:,n(1),j,k)+curl(:,n(1)+1,j,k)+curl(:,n(1),j+1,k)+curl(:,n(1)+1,j+1,k))/4
      endif
      ! Physical E=+eta*curl(B); original LFM stores minus physical E.
      emf(i,j,k,d)=emf(i,j,k,d)-w*vd*dl*dot_product(jedge,edge)
     enddo
    enddo
   enddo
  enddo
 end subroutine
 subroutine kaiju_nudge_outer(c,q,dt)
  class(cube_model),intent(in)::c
  real(real64),intent(inout)::q(1-ng:,1-ng:,1-ng:,:)
  real(real64),intent(in)::dt
  real(real64)::x(3),w,wind(8),alpha,rho,pressure
  integer::i,j,k,n(3)
  if(c%mode<2.or.c%layout%coord(1)/=c%layout%dims(1)-1) return
  n=c%n;i=n(1)
  do k=1,n(3)
   do j=1,n(2)
    x=c%outer_cell_xyz(:,0,j,k)
    w=c%outer_cell_weight(0,j,k)
    wind=earth_wind(c%sw,c%time,x,c%wind_switch_s,c%wind_bz_after,c%wind_reference_offset_s)
    alpha=min(1._real64,dt*w/60._real64)
    pressure=(1-alpha)*q(i,j,k,1)*q(i,j,k,5)**2/gamma+alpha*wind(1)*wind(5)**2/gamma
    rho=(1-alpha)*q(i,j,k,1)+alpha*wind(1)
    call demand(rho>0.and.pressure>0,'positive nudged boundary state')
    q(i,j,k,1)=rho;q(i,j,k,2:4)=(1-alpha)*q(i,j,k,2:4)+alpha*wind(2:4)
    q(i,j,k,5)=sqrt(gamma*pressure/rho)
   enddo
  enddo
 end subroutine
end module
