! Donor curl/project-before-interpolation on the ORIGINAL Cube mesh.
! Only geometry/index/ownership adapters are new; no new current or tracer.
module lfm_cube_yy_fac_source
 use iso_fortran_env,only:real64
 use iso_c_binding,only:c_int,c_double,c_ptr,c_null_ptr,c_associated
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 use mpi_f08
 use lfm_cube_mesh
 use lfm_cube_topology
 use lfm_cube_mpi,only:cube_halo_plan
 use lfm_cube_yy_fac_map,only:yy_released_source_map
 use lfm_cube_yy_igrf,only:cube_yy_igrf_field
 implicit none
 private
 public::yy_source_prepare,yy_source_refresh,yy_source_map_fac
 type(cube_meridian),save::mesh
 type(c_ptr),save::handle=c_null_ptr,ct_handle=c_null_ptr
 type(MPI_Comm),save::inner_comm=MPI_COMM_NULL
 integer,save::nlocal(3)=0
 integer,save::nang=0,njloc=0,nkloc=0,rank=-1,inner_count=0
 logical,save::inner=.false.,prepared=.false.,source_ready=.false.
 real(real64),save::length_m=0,source_epoch=-1
 integer,allocatable,save::descriptors(:,:)
 real(real64),allocatable,save::bpack(:,:,:,:),cart(:,:,:,:),parallel(:,:,:),tiles(:,:,:,:,:),moment_packet(:,:,:,:)
 real(real64),allocatable,save::global_parallel(:,:,:,:,:),ct_field(:,:,:,:),source_q(:,:,:,:)
 interface
  integer(c_int)function prepare_ct(ni,nj,nk,l,v,h)bind(C,name='cube_yy_ct_source_prepare')
   import c_int,c_double,c_ptr
   integer(c_int),value::ni,nj,nk
   real(c_double),value::l
   real(c_double),intent(in)::v(*)
   type(c_ptr),intent(out)::h
  end function
  integer(c_int)function apply_ct(h,p,bg,b)bind(C,name='cube_yy_ct_source_apply')
   import c_int,c_double,c_ptr
   type(c_ptr),value::h
   real(c_double),intent(in)::p(*),bg(*)
   real(c_double),intent(out)::b(*)
  end function
  integer(c_int)function prepare_current(nj,nk,l,v,h)bind(C,name='cube_yy_current_prepare')
   import c_int,c_double,c_ptr
   integer(c_int),value::nj,nk
   real(c_double),value::l
   real(c_double),intent(in)::v(*)
   type(c_ptr),intent(out)::h
  end function
  integer(c_int)function apply_current(h,b,t,j,p)bind(C,name='cube_yy_current_apply')
   import c_int,c_double,c_ptr
   type(c_ptr),value::h
   real(c_double),intent(in)::b(*)
   real(c_double),value::t
   real(c_double),intent(out)::j(*),p(*)
  end function
 end interface
contains
 subroutine collective_status(comm,local,status)
  type(MPI_Comm),intent(in)::comm
  integer,intent(in)::local
  integer,intent(out)::status
  integer::st
  call MPI_Allreduce(local,status,1,MPI_INTEGER,MPI_MAX,comm,st)
  if(st/=MPI_SUCCESS)status=99
 end subroutine
 subroutine yy_source_prepare(m,layout,g,vertices,model_length_m,comm,status)
  type(cube_meridian),intent(in)::m
  type(cube_layout),intent(in)::layout
  integer,intent(in)::g
  real(real64),intent(in)::vertices(:,-g:,-g:,-g:),model_length_m
  type(MPI_Comm),intent(in)::comm
  integer,intent(out)::status
  real(real64),allocatable::v(:,:,:,:),vct(:,:,:,:)
  integer,allocatable::coverage(:,:,:)
  integer::colour,st,local,descriptor(4),p,j,k,lo(2),hi(2),f,j0,k0
  ! One model per MPI process; caches are rebuilt by normal own-restart init.
  local=0;if(prepared.or.g<3.or.model_length_m<=0) local=1
  if(layout%global(2)/=layout%global(3).or.layout%local(1)<6) local=1
  call collective_status(comm,local,status);if(status/=0)return
  mesh=m;length_m=model_length_m;nang=layout%global(2);nlocal=layout%local
  ! Exact generic donor CT recovery; each MPI owner prepares its own active grid once.
  vct=vertices(:,0:nlocal(1),0:nlocal(2),0:nlocal(3))
  local=prepare_ct(int(nlocal(1),c_int),int(nlocal(2),c_int),int(nlocal(3),c_int),length_m,vct,ct_handle)
  call collective_status(comm,local,status);if(status/=0)return
  allocate(ct_field(nlocal(1),nlocal(2),nlocal(3),3))
  allocate(source_q(1-g:nlocal(1)+g,1-g:nlocal(2)+g,1-g:nlocal(3)+g,8))
  njloc=layout%local(2);nkloc=layout%local(3);rank=layout%rank
  inner=layout%coord(1)==0
  colour=MPI_UNDEFINED;if(inner)colour=1
  call MPI_Comm_split(comm,colour,rank,inner_comm,st)
  local=0;if(st/=MPI_SUCCESS)local=99
  call collective_status(comm,local,status);if(status/=0)return
  allocate(global_parallel(0:5,0:nang+1,0:nang+1,0:5,3));global_parallel=0
  if(inner)then
   v=vertices(:,-2:6,-2:njloc+2,-2:nkloc+2)
   local=prepare_current(int(njloc,c_int),int(nkloc,c_int),length_m,v,handle)
   if(local==0)then
    allocate(bpack(3,-1:6,-1:njloc+2,-1:nkloc+2),cart(3,0:5,0:njloc+1,0:nkloc+1))
    allocate(parallel(0:5,0:njloc+1,0:nkloc+1))
    allocate(moment_packet(0:5,0:njloc+1,0:nkloc+1,3))
   endif
   call MPI_Comm_size(inner_comm,inner_count,st)
   if(st/=MPI_SUCCESS)local=99
  endif
  call collective_status(comm,local,status);if(status/=0)return
  if(inner)then
   if(rank==0)then
    allocate(descriptors(4,inner_count),tiles(0:5,0:njloc+1,0:nkloc+1,3,inner_count))
   else
    allocate(descriptors(4,1),tiles(1,1,1,1,1))
   endif
   descriptor=[layout%face,layout%first(2),layout%first(3),rank]
   call MPI_Gather(descriptor,4,MPI_INTEGER,descriptors,4,MPI_INTEGER,0,inner_comm,st)
   if(st/=MPI_SUCCESS)local=99
  endif
  call collective_status(comm,local,status);if(status/=0)return
  if(rank==0)then
   allocate(coverage(0:nang+1,0:nang+1,0:5));coverage=0
   do p=1,inner_count
    f=descriptors(1,p);j0=descriptors(2,p);k0=descriptors(3,p)
    lo=[1,1];hi=[njloc,nkloc]
    if(j0==0)lo(1)=0
    if(k0==0)lo(2)=0
    if(j0+njloc==nang)hi(1)=njloc+1
    if(k0+nkloc==nang)hi(2)=nkloc+1
    do k=lo(2),hi(2)
     do j=lo(1),hi(1)
      coverage(j0+j,k0+k,f)=coverage(j0+j,k0+k,f)+1
     enddo
    enddo
   enddo
   if(any(coverage/=1))local=2
  endif
  call collective_status(comm,local,status);if(status/=0)return
  prepared=.true.;source_ready=.false.
 end subroutine
 subroutine yy_source_refresh(q,phi,bgphi,halo,g,time,comm,status)
  integer,intent(in)::g
  real(real64),intent(in)::q(1-g:,1-g:,1-g:,:),time
  real(real64),contiguous,intent(in)::phi(0:,0:,0:,:),bgphi(0:,0:,0:,:)
  type(cube_halo_plan),intent(inout)::halo
  logical::bounded(8),positive(8)
  type(MPI_Comm),intent(in)::comm
  integer,intent(out)::status
  integer::i,j,k,p,f,j0,k0,lo(2),hi(2),st,local
  local=0;if(.not.prepared.or..not.ieee_is_finite(time))local=1
  call collective_status(comm,local,status);if(status/=0)return
  ! Recover from actual independent CT flux. Archive total-cellB is never used.
  local=apply_ct(ct_handle,phi,bgphi,ct_field)
  call collective_status(comm,local,status);if(status/=0)return
  source_q=q
  do k=1,nlocal(3)
   do j=1,nlocal(2)
    do i=1,nlocal(1)
     source_q(i,j,k,6:8)=ct_field(i,j,k,:)
    enddo
   enddo
  enddo
  bounded=.false.;positive=.false.;positive([1,5])=.true.
  call halo%sync(source_q,st,bounded,positive)
  local=0;if(st/=0)local=99
  call collective_status(comm,local,status);if(status/=0)return
  if(inner)then
   ! Same strict donor cell-B ghost rule. The source layer is cell 2, so its
   ! curl uses active neighbours, rather than this first-cell wall stencil.
   do i=1-g,0
    source_q(i,:,:,:)=source_q(1,:,:,:)
   enddo
   do k=-1,nkloc+2
    do j=-1,njloc+2
     do i=-1,6
      bpack(1,i,j,k)=source_q(i,j,k,6)
      bpack(2,i,j,k)=source_q(i,j,k,7)
      bpack(3,i,j,k)=source_q(i,j,k,8)
     enddo
    enddo
   enddo
   local=apply_current(handle,bpack,time,cart,parallel)
  endif
  call collective_status(comm,local,status);if(status/=0)return
  if(inner)then
   moment_packet(:,:,:,1)=parallel
   moment_packet(:,:,:,2)=source_q(0:5,0:njloc+1,0:nkloc+1,1)
   moment_packet(:,:,:,3)=source_q(0:5,0:njloc+1,0:nkloc+1,5)
   call MPI_Gather(moment_packet,size(moment_packet),MPI_DOUBLE_PRECISION,tiles,size(moment_packet), &
    MPI_DOUBLE_PRECISION,0,inner_comm,st)
   if(st/=MPI_SUCCESS)local=99
  endif
  call collective_status(comm,local,status);if(status/=0)return
  if(rank==0)then
   do p=1,inner_count
    f=descriptors(1,p);j0=descriptors(2,p);k0=descriptors(3,p)
    lo=[1,1];hi=[njloc,nkloc]
    if(j0==0)lo(1)=0
    if(k0==0)lo(2)=0
    if(j0+njloc==nang)hi(1)=njloc+1
    if(k0+nkloc==nang)hi(2)=nkloc+1
    do k=lo(2),hi(2)
     do j=lo(1),hi(1)
      global_parallel(:,j0+j,k0+k,f,:)=tiles(:,j,k,:,p)
     enddo
    enddo
   enddo
  endif
  call MPI_Bcast(global_parallel,size(global_parallel),MPI_DOUBLE_PRECISION,0,comm,st)
  local=0;if(st/=MPI_SUCCESS)local=99
  if(.not.all(ieee_is_finite(global_parallel)))local=3
  call collective_status(comm,local,status);if(status/=0)return
  source_epoch=time;source_ready=.true.
 end subroutine
 subroutine yy_source_map_fac(x,y,theta,radii,rion_m,time,rn,rs,trace,status,work_first,work_last,moments)
  real(real64),intent(in)::x(:,:),y(:,:),theta(:),radii(:),rion_m,time
  real(real64),intent(out)::rn(:,:),rs(:,:)
  integer,intent(out)::trace(:,:,:),status
  integer,intent(in)::work_first,work_last
  real(real64),intent(out)::moments(:,:,:,:)
  rn=0;rs=0;trace=-999;status=1
  if(.not.source_ready.or.time/=source_epoch)return
  call yy_released_source_map(mesh,length_m,global_parallel,x,y,theta,radii,rion_m,time, &
   rn,rs,trace,status,work_first,work_last,moments)
 end subroutine
end module
