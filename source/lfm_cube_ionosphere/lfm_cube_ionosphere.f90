! Cube thin-shell transport around the unchanged original LFM ionosphere.
! No elliptic operator, FAC curl, conductance model or field-line solver here.
module lfm_cube_ionosphere
 use iso_fortran_env,only:real64,int64
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 use mpi_f08
 use lfm_cube_topology
 use lfm_cube_mesh
 use lfm_cube_mpi,only:polynomial_weights
 use lfm_cube_yy_fac_source,only:yy_source_map_fac
 use lfm_cube_yy_boundary,only:yy_boundary_install,yy_boundary_potential, &
  yy_boundary_collective_refresh,yy_boundary_volume_domain,yy_boundary_set_time
 implicit none
 private
 public::cube_ionosphere
 public::cube_mi_nj,cube_mi_nk,cube_mi_capj,cube_mi_capk,cube_mi_state_count
 include 'native_mi_shape.inc'
 include 'yy_rotation_profile.inc'
 integer,parameter::cube_mi_nj=nj,cube_mi_nk=nk,cube_mi_capj=capj,cube_mi_capk=capk
 integer,parameter::cube_mi_state_count=header+22*capj*capk+capj*(nk+1)
 real(real64),parameter::rion=6.5e6_real64
 type cube_ionosphere
  type(cube_layout)::layout
  type(MPI_Comm)::comm=MPI_COMM_NULL,inner_comm=MPI_COMM_NULL
  integer::rank=-1,nang=0,last_step=-1,nunique=0
  real(real64)::last_time=-1,sigma=5,worst_residual=0,cpcp(2)=0,facmax=0
  real(real64)::gather_seconds=0,sample_seconds=0,solve_seconds=0,wall_seconds=0
  real(real64)::maximum_moment_residual=0,maximum_weight_norm=0,maximum_duplicate_point_error=0
  logical::split_background=.false.,native_primed=.false.,ready=.false.
  integer,allocatable::donor(:,:,:),vertex(:,:,:),inner_tiles(:,:)
  real(real64),allocatable::weight(:,:),reference_xyz(:,:,:,:),reference_background(:,:,:,:)
  real(real64),allocatable::theta(:),azimuth(:),wall_xyz(:,:),wall_phi(:)
  real(real64)::pn(capj,capk)=0,ps(capj,capk)=0,table(capj,nk+1)=0
  ! Exact native snapshot: ten documented fields, north then south.
  real(real64)::spatial(capj,capk,10,2)=0
  real(real64)::cap_x(capj,capk)=0,cap_y(capj,capk)=0,source_radii(capj)=0
  real(real64)::source_theta(capj)=0,source_phi(nk+1)=0,solve_SER(2)=0
  integer::solve_ITER(2)=0,solve_converged(2)=1,solve_ITMAX(2)=100
  integer::mapped_trace_status(capj,capk,2)=-999
 contains
  procedure::init=>ion_init
  procedure::sample=>ion_sample
  procedure::update=>ion_update
  procedure::wall_values=>ion_wall_values
  procedure::snapshot=>ion_snapshot
  procedure::export=>ion_export
  procedure::import=>ion_import
  procedure::checkpoint_export=>ion_checkpoint_export
  procedure::checkpoint_import=>ion_checkpoint_import
 end type
 interface
  subroutine lfm_cube_mi_setup(geometry,count,pedersen,status)
   import real64
   integer,intent(in)::count
   real(real64),intent(in)::geometry(count),pedersen
   integer,intent(out)::status
  end subroutine
  subroutine lfm_sw_ionosphere_update_packed(step,time,dt,state,count,status,north,south,fac)
   import real64
   integer,intent(in)::step,count
   real(real64),intent(in)::time,dt,state(count)
   integer,intent(out)::status
   real(real64),intent(out)::north,south,fac
  end subroutine
  subroutine lfm_cube_mi_get(pn,ps,table,worst,status)
   import real64,capj,capk,nk
   real(real64),intent(out)::pn(capj,capk),ps(capj,capk),table(capj,nk+1),worst
   integer,intent(out)::status
  end subroutine
  subroutine lfm_cube_mi_put(pn,ps,status)
   import real64,capj,capk
   real(real64),intent(in)::pn(capj,capk),ps(capj,capk)
   integer,intent(out)::status
  end subroutine
  subroutine lfm_cube_mi_restore_held(pn,ps,cn,cs,step,time,table,spatial,status)
   import real64,capj,capk,nk
   real(real64),intent(in)::pn(capj,capk),ps(capj,capk),cn(capj,capk),cs(capj,capk),time
   integer,intent(in)::step
   real(real64),intent(out)::table(capj,nk+1),spatial(capj,capk,10,2)
   integer,intent(out)::status
  end subroutine
  subroutine lfm_sw_mi_snapshot_pack(nj_out,nk_out,data,status)
   import real64
   integer,intent(in)::nj_out,nk_out
   real(real64),intent(out)::data(nj_out,nk_out,10,2)
   integer,intent(out)::status
  end subroutine
 end interface
contains
 subroutine collective_status(comm,local,status)
  type(MPI_Comm),intent(in)::comm
  integer,intent(in)::local
  integer,intent(out)::status
  integer::ierr
  call MPI_Allreduce(local,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  if(ierr/=MPI_SUCCESS) status=99
 end subroutine
 subroutine ion_init(h,m,layout,reference_grid,comm,status,sigma,ref_background_B)
  class(cube_ionosphere),intent(out)::h
  type(cube_meridian),intent(in)::m
  type(cube_layout),intent(in)::layout
  character(*),intent(in)::reference_grid
  type(MPI_Comm),intent(in)::comm
  integer,intent(out)::status
  real(real64),optional,intent(in)::sigma,ref_background_B(:,:,:,:)
  character(80)::magic
  integer::local,st,u,dims(3),i,j,k,d,n,colour,split,ntiles,t,descriptor(5),f,j0,k0,njloc,nkloc,comm_ranks
  integer,allocatable::coverage(:,:,:)
  real(real64),allocatable::packed(:)
  h%comm=comm;h%layout=layout;h%nang=layout%global(2)
  call MPI_Comm_rank(comm,h%rank,st)
  call MPI_Comm_size(comm,comm_ranks,st)
  if(present(sigma)) h%sigma=sigma
  local=0
  if(any(layout%dims<1).or.h%nang<4.or.layout%local(1)<nc) local=1
  if(6*product(layout%dims)/=layout%size.or.any(mod(layout%global,max(1,layout%dims))/=0)) local=1
  if(layout%rank/=h%rank.or.layout%size/=comm_ranks) local=1
  if(layout%global(2)/=layout%global(3).or.layout%dims(2)/=layout%dims(3)) local=1
  ! The reviewed native/Cube pairs use matching angular refinement. Reject
  ! a fine MHD shell linked to the coarse native MI wall (or the reverse).
  if(nj/=2*h%nang.or.3*nk/=8*h%nang) local=1
  if(h%sigma<=0.or..not.ieee_is_finite(h%sigma)) local=1
  call collective_status(comm,local,status);if(status/=0) return
  allocate(h%reference_xyz(3,5,capj,nk+1))
  if(h%rank==0) then
   open(newunit=u,file=reference_grid,status='old',action='read',iostat=local)
   if(local==0) then
    read(u,*,iostat=local) magic
    if(local==0) read(u,*,iostat=local) dims
    if(local==0) then
     if(magic/='LFM_CUBE_ION_REFERENCE_V1'.or.any(dims/=[5,capj,nk+1])) local=1
    endif
    do k=1,nk+1
     do j=1,capj
      do i=1,5
       if(local/=0) exit
       read(u,*,iostat=local) h%reference_xyz(:,i,j,k)
      enddo
     enddo
    enddo
    close(u)
   endif
   if(local==0) then
    if(.not.all(ieee_is_finite(h%reference_xyz))) local=1
   endif
  endif
  call MPI_Bcast(local,1,MPI_INTEGER,0,comm,st);status=abs(local);if(status/=0) return
  call MPI_Bcast(h%reference_xyz,size(h%reference_xyz),MPI_DOUBLE_PRECISION,0,comm,st)
  split=0;if(present(ref_background_B)) split=1
  call MPI_Bcast(split,1,MPI_INTEGER,0,comm,st)
  local=0
  if((split==1).neqv.present(ref_background_B)) local=1
  if(present(ref_background_B)) then
   if(any(shape(ref_background_B)/=[4,nj,nk,3])) local=1
   if(.not.all(ieee_is_finite(ref_background_B))) local=1
   if(h%rank==0.and.local==0) then
    allocate(h%reference_background(4,nj,nk,3));h%reference_background=ref_background_B
   endif
  endif
  h%split_background=split==1
  call collective_status(comm,local,status);if(status/=0) return
  colour=MPI_UNDEFINED;if(layout%coord(1)==0) colour=1
  call MPI_Comm_split(comm,colour,layout%rank,h%inner_comm,st)
  ! Only wall-owning tiles participate. Reassemble each global face from its
  ! angular tiles without gathering any interior radial state.
  ntiles=6*layout%dims(2)*layout%dims(3)
  n=1;if(h%rank==0) n=ntiles
  allocate(h%inner_tiles(5,n))
  descriptor=[layout%face,layout%first(2),layout%first(3),layout%local(2),layout%local(3)]
  if(layout%coord(1)==0) call MPI_Gather(descriptor,5,MPI_INTEGER,h%inner_tiles,5,MPI_INTEGER,0,h%inner_comm,st)
  local=0
  if(h%rank==0) then
   allocate(coverage(h%nang,h%nang,0:5));coverage=0
   do t=1,ntiles
    f=h%inner_tiles(1,t);j0=h%inner_tiles(2,t);k0=h%inner_tiles(3,t)
    njloc=h%inner_tiles(4,t);nkloc=h%inner_tiles(5,t)
    if(f<0.or.f>5.or.min(j0,k0)<0.or.min(njloc,nkloc)<=0.or. &
       max(j0+njloc,k0+nkloc)>h%nang) then
     local=1;cycle
    endif
    coverage(j0+1:j0+njloc,k0+1:k0+nkloc,f)=coverage(j0+1:j0+njloc,k0+1:k0+nkloc,f)+1
   enddo
   if(any(coverage/=1)) local=1
  endif
  call collective_status(comm,local,status);if(status/=0) return
  allocate(h%theta(capj),h%azimuth(nk+1))
  h%theta=atan2(sqrt(h%reference_xyz(2,1,:,1)**2+h%reference_xyz(3,1,:,1)**2),h%reference_xyz(1,1,:,1))
  do k=1,nk+1
   h%azimuth(k)=atan2(h%reference_xyz(3,1,2,k),h%reference_xyz(2,1,2,k))
   if(k==1) cycle
   do while(h%azimuth(k)<=h%azimuth(k-1))
    h%azimuth(k)=h%azimuth(k)+2*acos(-1._real64)
   enddo
  enddo
  local=0
  if(h%rank==0) call build_sampler(h,m,local)
  call MPI_Bcast(local,1,MPI_INTEGER,0,comm,st);status=abs(local);if(status/=0) return
  call build_wall_vertices(h,m,local)
  call collective_status(comm,local,status);if(status/=0) return
  call configure_return_domain(h,m,status);if(status/=0) return
  if(h%rank==0) then
   allocate(packed(3*4*capj*(nk+1)));n=0
   do d=1,3
    packed(n+1:n+4*capj*(nk+1))=reshape(h%reference_xyz(d,1:4,:,:),[4*capj*(nk+1)])
    n=n+4*capj*(nk+1)
   enddo
   call lfm_cube_mi_setup(packed,size(packed),h%sigma,local)
   if(local==0) call lfm_cube_yy_mi100_policy(local)
   if(local==0) call lfm_cube_yy_cap_geometry(h%cap_x,h%cap_y,h%source_theta,h%source_radii,h%source_phi,local)
  endif
  call MPI_Bcast(local,1,MPI_INTEGER,0,comm,st);status=abs(local)
  call MPI_Bcast(h%cap_x,size(h%cap_x),MPI_DOUBLE_PRECISION,0,comm,st)
  call MPI_Bcast(h%cap_y,size(h%cap_y),MPI_DOUBLE_PRECISION,0,comm,st)
  if(st/=MPI_SUCCESS) status=99
  call MPI_Bcast(h%source_theta,size(h%source_theta),MPI_DOUBLE_PRECISION,0,comm,st)
  if(st/=MPI_SUCCESS) status=99
  call MPI_Bcast(h%source_phi,size(h%source_phi),MPI_DOUBLE_PRECISION,0,comm,st)
  if(st/=MPI_SUCCESS) status=99
  call MPI_Bcast(h%source_radii,size(h%source_radii),MPI_DOUBLE_PRECISION,0,comm,st)
  if(st/=MPI_SUCCESS) status=99
  h%ready=status==0
 end subroutine
 subroutine build_sampler(h,m,status)
  class(cube_ionosphere),intent(inout)::h
  type(cube_meridian),intent(in)::m
  integer,intent(out)::status
  real(real64)::v(3,0:1,0:1,0:1),rv(3,0:1,0:1,0:1),centre(3),cov(3,3),volume
  real(real64)::targetcentre(3),targetcov(3,3),a,b,t,direction(3),delta(3,3),transform(3,3)
  real(real64)::origin(3),matrix(27,10),target(10),w(27),residual
  integer::i,j,k,q,f,idx(3),start(3),ii,jj,kk,z,st,c
  status=1;allocate(h%donor(4,27,nc*nj*nk),h%weight(27,nc*nj*nk));q=0
  do k=1,nk
   do j=1,nj
    do i=1,nc
     q=q+1
     do c=0,7
      rv(:,mod(c,2),mod(c/2,2),c/4)=h%reference_xyz(:,i+mod(c,2),j+mod(c/2,2),k+c/4)
     enddo
     call cube_cell_statistics(rv,targetcentre,targetcov,volume,st);if(st/=0) return
     call m%reference(targetcentre,t,direction,st);if(st/=0) return
     call cube_inverse(direction,f,a,b,st);if(st/=0) return
     idx=[int(floor(t*m%ni))+1,int(floor((a/(acos(-1._real64)/4)+1)*h%nang/2))+1, &
      int(floor((b/(acos(-1._real64)/4)+1)*h%nang/2))+1]
     idx=max(1,min([nc,h%nang,h%nang],idx));start=max(1,min([nc-2,h%nang-2,h%nang-2],idx-1))
     call m%cell_vertices(f,h%nang,idx,v,st);if(st/=0) return
     delta(:,1)=sum(sum(v(:,1,:,:)-v(:,0,:,:),dim=3),dim=2)/4
     delta(:,2)=sum(sum(v(:,:,1,:)-v(:,:,0,:),dim=3),dim=2)/4
     delta(:,3)=sum(sum(v(:,:,:,1)-v(:,:,:,0),dim=3),dim=2)/4
     call cube_inverse3(delta,transform,st);if(st/=0) return
     origin=sum(sum(sum(v,dim=4),dim=3),dim=2)/8;z=0
     do kk=0,2
      do jj=0,2
       do ii=0,2
        z=z+1;h%donor(:,z,q)=[f,start+[ii,jj,kk]]
        call m%cell_vertices(f,h%nang,start+[ii,jj,kk],v,st);if(st/=0) return
        call cube_cell_statistics(v,centre,cov,volume,st);if(st/=0) return
        matrix(z,:)=cube_polynomial_moments(centre,cov,origin,transform)
       enddo
      enddo
     enddo
     target=cube_polynomial_moments(targetcentre,targetcov,origin,transform)
     call polynomial_weights(matrix,target,w,residual,st);if(st/=0) return
     h%weight(:,q)=w
     h%maximum_moment_residual=max(h%maximum_moment_residual,residual)
     h%maximum_weight_norm=max(h%maximum_weight_norm,sum(abs(w)))
    enddo
   enddo
  enddo
  status=0
 end subroutine
 subroutine configure_return_domain(h,m,status)
  ! Query the existing mapped geometry at the very same inner radial index
  ! and angular ghost range used by model_guards (ng=4). No geometry fit.
  class(cube_ionosphere),intent(in)::h
  type(cube_meridian),intent(in)::m
  integer,intent(out)::status
  integer,parameter::wall_ng=4
  real(real64),allocatable::xyz(:,:,:)
  real(real64)::bounds(2),global_bounds(2),x(3),r
  integer::j,k,st,local,ierr,njloc,nkloc
  bounds=[-huge(1._real64),0._real64];local=0
  if(h%layout%coord(1)==0) then
   njloc=h%layout%local(2);nkloc=h%layout%local(3)
   allocate(xyz(3,-wall_ng:njloc+wall_ng,-wall_ng:nkloc+wall_ng))
   do k=-wall_ng,nkloc+wall_ng
    do j=-wall_ng,njloc+wall_ng
     call m%point(h%layout%face,h%nang,h%layout%first+[0,j,k],x,st)
     if(st/=0.or..not.all(ieee_is_finite(x))) then
      local=1;cycle
     endif
     xyz(:,j,k)=x;r=norm2(x);bounds=max(bounds,[-r,r])
    enddo
   enddo
   if(local==0) then
    do k=1-wall_ng,nkloc+wall_ng
     do j=1-wall_ng,njloc+wall_ng
      x=(xyz(:,j-1,k-1)+xyz(:,j,k-1)+xyz(:,j,k)+xyz(:,j-1,k))/4
      r=norm2(x);bounds=max(bounds,[-r,r])
     enddo
    enddo
   endif
  endif
  call collective_status(h%comm,local,status);if(status/=0) return
  call MPI_Allreduce(bounds,global_bounds,2,MPI_DOUBLE_PRECISION,MPI_MAX,h%comm,ierr)
  if(ierr/=MPI_SUCCESS) then
   status=99;return
  endif
  call yy_boundary_volume_domain(-global_bounds(1),global_bounds(2),rion,local)
  call collective_status(h%comm,local,status)
 end subroutine
 subroutine build_wall_vertices(h,m,status)
  class(cube_ionosphere),intent(inout)::h
  type(cube_meridian),intent(in)::m
  integer,intent(out)::status
  integer(int64),allocatable::keys(:,:)
  integer(int64)::key(4)
  integer::f,j,k,l,found,st,capacity
  real(real64)::x(3),error
  status=1;capacity=6*(h%nang+1)**2
  allocate(h%vertex(0:h%nang,0:h%nang,0:5),keys(4,capacity),h%wall_xyz(3,capacity))
  h%nunique=0
  do f=0,5
   do k=0,h%nang
    do j=0,h%nang
     call cube_vertex_key(f,h%layout%global,[0,j,k],key,st);if(st/=0) return
     call m%point(f,h%nang,[0,j,k],x,st);if(st/=0) return
     found=0
     do l=1,h%nunique
      if(all(keys(:,l)==key)) then
       found=l;exit
      endif
     enddo
     if(found==0) then
      h%nunique=h%nunique+1;found=h%nunique;keys(:,found)=key;h%wall_xyz(:,found)=x
     else
      error=norm2(x-h%wall_xyz(:,found))/max(1._real64,norm2(x))
      h%maximum_duplicate_point_error=max(h%maximum_duplicate_point_error,error)
      if(error>1.e-10_real64) return
     endif
     h%vertex(j,k,f)=found
    enddo
   enddo
  enddo
  if(h%nunique/=6*h%nang**2+2) return
  allocate(h%wall_phi(h%nunique));h%wall_phi=0;status=0
 end subroutine
 subroutine ion_sample(h,inner,reference,status)
  class(cube_ionosphere),intent(inout)::h
  real(real64),intent(in)::inner(:,:,:,:)
  real(real64),allocatable,intent(out)::reference(:,:,:,:)
  integer,intent(out)::status
  real(real64),allocatable::received(:,:,:,:,:),tiles(:,:,:,:,:)
  real(real64)::started,value,lo
  integer::local,st,q,z,i,j,k,c,index(4),ntiles,t,f,j0,k0,njloc,nkloc
  status=1;if(.not.h%ready) return
  local=0
  if(h%layout%coord(1)==0) then
   if(any(shape(inner)/=[nc,h%layout%local(2),h%layout%local(3),8])) local=1
   if(.not.all(ieee_is_finite(inner))) local=1
   if(size(inner,4)==8) then
    if(minval(inner(:,:,:,1))<=0.or.minval(inner(:,:,:,5))<=0) local=1
   endif
  endif
  call collective_status(h%comm,local,status);if(status/=0) return
  started=MPI_Wtime()
  ntiles=6*h%layout%dims(2)*h%layout%dims(3)
  if(h%rank==0) then
   allocate(received(nc,h%nang,h%nang,8,6),reference(nc,nj,nk,8))
   allocate(tiles(nc,h%layout%local(2),h%layout%local(3),8,ntiles))
  else
   allocate(received(1,1,1,1,1),reference(0,0,0,0))
   allocate(tiles(1,1,1,1,1))
  endif
  if(h%layout%coord(1)==0) then
   call MPI_Gather(inner,size(inner),MPI_DOUBLE_PRECISION,tiles,size(inner),MPI_DOUBLE_PRECISION,0,h%inner_comm,st)
  endif
  h%gather_seconds=h%gather_seconds+MPI_Wtime()-started;started=MPI_Wtime();local=0
  if(h%rank==0) then
   do t=1,ntiles
    f=h%inner_tiles(1,t);j0=h%inner_tiles(2,t);k0=h%inner_tiles(3,t)
    njloc=h%inner_tiles(4,t);nkloc=h%inner_tiles(5,t)
    received(:,j0+1:j0+njloc,k0+1:k0+nkloc,:,f+1)=tiles(:,:,:,:,t)
   enddo
   q=0
   do k=1,nk
    do j=1,nj
     do i=1,nc
      q=q+1
      do c=1,8
       value=0
       do z=1,27
        index=h%donor(:,z,q)
        value=value+h%weight(z,q)*received(index(2),index(3),index(4),c,index(1)+1)
       enddo
        if(c==1.or.c==5) then
         ! Same positive-scalar transfer contract as cube_halo_plan%sync:
         ! reject invalid donors; a signed high-order undershoot uses their
         ! positive lower bound. This changes only the MI reference packet,
         ! never the live MHD density, sound speed, CFL or physical floors.
         lo=huge(value)
         do z=1,27
          index=h%donor(:,z,q)
          lo=min(lo,received(index(2),index(3),index(4),c,index(1)+1))
         enddo
         if(lo<=0) local=1
         if(value<=0) value=lo
        endif
       reference(i,j,k,c)=value
      enddo
     enddo
    enddo
   enddo
   if(h%split_background) reference(:,:,:,6:8)=reference(:,:,:,6:8)+h%reference_background
   if(.not.all(ieee_is_finite(reference))) local=1
   if(minval(reference(:,:,:,1))<=0.or.minval(reference(:,:,:,5))<=0) local=1
  endif
  h%sample_seconds=h%sample_seconds+MPI_Wtime()-started
  call MPI_Bcast(local,1,MPI_INTEGER,0,h%comm,st);status=local
 end subroutine
 subroutine evaluate_wall(h,status)
  class(cube_ionosphere),intent(inout)::h
  integer,intent(out)::status
  integer::i,st,local,first,last,p,ranks
  integer::counts(h%layout%size),offsets(h%layout%size)
  ranks=h%layout%size;local=0
  do p=0,ranks-1
   offsets(p+1)=h%nunique*p/ranks
   counts(p+1)=h%nunique*(p+1)/ranks-offsets(p+1)
  enddo
  first=offsets(h%rank+1)+1;last=first+counts(h%rank+1)-1
  do i=first,last
   call yy_boundary_potential(.01_real64*h%wall_xyz(:,i),h%wall_phi(i),st)
   if(st/=0) then
    write(9,*) 'YY_MI_WALL_POINT_FAIL',i,st,h%wall_xyz(:,i);flush(9)
    local=1;exit
   endif
  enddo
  call collective_status(h%comm,local,status);if(status/=0) return
  call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DOUBLE_PRECISION,h%wall_phi,counts,offsets, &
   MPI_DOUBLE_PRECISION,h%comm,st)
  if(st/=MPI_SUCCESS) status=99
 end subroutine
 subroutine ion_update(h,inner,step,time,status,cpcp,facmax,dt_s)
  class(cube_ionosphere),intent(inout)::h
  real(real64),intent(in)::inner(:,:,:,:),time
  integer,intent(in)::step
  integer,intent(out)::status
  real(real64),optional,intent(out)::cpcp(2),facmax
  real(real64),optional,intent(in)::dt_s
  real(real64),allocatable::reference(:,:,:,:),packed(:)
  real(real64)::started,dt,jcart(3,capj,nk+1),rn(capj,capk),rs(capj,capk),unused(3)
  integer::st,d,n,nfield,local,capture_unit,capture_io
  integer::p,first,last,total,ranks
  integer::counts(h%layout%size),offsets(h%layout%size),trace_flat(2*capj*capk),trace_all(2*capj*capk)
  real(real64)::mapped_flat(2*capj*capk),mapped_all(2*capj*capk)
  local=0
  if(step<0.or.time<0.or..not.ieee_is_finite(time)) local=1
  if(step<h%last_step.or.time<h%last_time) local=1
  dt=0
  if(present(dt_s)) dt=dt_s
  if(dt<0.or..not.ieee_is_finite(dt)) local=1
  call collective_status(h%comm,local,status);if(status/=0) return
  call h%sample(inner,reference,status);if(status/=0) return
  local=0;started=MPI_Wtime()
  if(h%rank==0) then
   allocate(packed((2*3+3*4)*nj*nk));n=0
   do d=1,5
    nfield=3;if(d>=3) nfield=4
    select case(d)
    case(1);packed(n+1:n+nfield*nj*nk)=reshape(reference(1:nfield,:,:,1),[nfield*nj*nk])
    case(2);packed(n+1:n+nfield*nj*nk)=reshape(reference(1:nfield,:,:,5),[nfield*nj*nk])
    case default;packed(n+1:n+nfield*nj*nk)=reshape(reference(1:nfield,:,:,d+3),[nfield*nj*nk])
    end select
    n=n+nfield*nj*nk
   enddo
   if(step==0) then
    open(newunit=capture_unit,file='yy_initial_mi_input.bin',access='stream',form='unformatted', &
     status='new',iostat=capture_io)
    if(capture_io==0) then
     write(capture_unit) 'YYMIP001',step,time,size(packed),packed,h%cap_x,h%cap_y, &
      h%source_theta,h%source_phi,h%source_radii
     close(capture_unit)
    endif
   endif
   call lfm_cube_yy_load_packed(step,time,dt,packed,size(packed),local,unused(1),unused(2),unused(3))
   write(9,*) 'YY_MI_STAGE LOAD',step,time,local;flush(9)
   write(9,*) 'YY_MI_STAGE RELEASED_SOURCE_READY',step,time,local;flush(9)
  endif
  ! Distribute independent cap points, preserving the original scalar trace,
  ! interpolation and multiplication order. Gatherv copies owner words; no
  ! floating-point sum or changed geometric approximation is introduced.
  call MPI_Bcast(local,1,MPI_INTEGER,0,h%comm,st)
  if(st/=MPI_SUCCESS) local=99
  call collective_status(h%comm,local,status);if(status/=0) return
  ranks=h%layout%size;total=2*capj*capk
  do p=0,ranks-1
   offsets(p+1)=total*p/ranks;counts(p+1)=total*(p+1)/ranks-offsets(p+1)
  enddo
  first=offsets(h%rank+1)+1;last=first+counts(h%rank+1)-1
  call yy_source_map_fac(h%cap_x,h%cap_y,h%source_theta,h%source_radii,rion,time, &
   rn,rs,h%mapped_trace_status,local,first,last)
  call collective_status(h%comm,local,status);if(status/=0) return
  mapped_flat=[reshape(rn,[capj*capk]),reshape(rs,[capj*capk])]
  trace_flat=reshape(h%mapped_trace_status,[total])
  call MPI_Gatherv(mapped_flat(first:last),counts(h%rank+1),MPI_DOUBLE_PRECISION, &
   mapped_all,counts,offsets,MPI_DOUBLE_PRECISION,0,h%comm,st)
  local=0;if(st/=MPI_SUCCESS) local=99
  call collective_status(h%comm,local,status);if(status/=0) return
  call MPI_Gatherv(trace_flat(first:last),counts(h%rank+1),MPI_INTEGER, &
   trace_all,counts,offsets,MPI_INTEGER,0,h%comm,st)
  local=0;if(st/=MPI_SUCCESS) local=99
  call collective_status(h%comm,local,status);if(status/=0) return
  if(h%rank==0) then
   rn=reshape(mapped_all(:capj*capk),shape(rn));rs=reshape(mapped_all(capj*capk+1:),shape(rs))
   h%mapped_trace_status=reshape(trace_all,shape(h%mapped_trace_status))
   write(9,*) 'YY_MI_STAGE MAP',step,time,local;flush(9)
   if(local==0) call lfm_cube_yy_solve_cap(rn,rs,h%pn,h%ps,local)
   write(9,*) 'YY_MI_STAGE SOLVE',step,time,local;flush(9)
   if(local==0) then
    call lfm_cube_yy_solve_metrics(h%solve_SER,h%solve_ITER,h%solve_ITMAX,h%solve_converged)
    if(all(rn==0)) then
     h%solve_SER(1)=0;h%solve_ITER(1)=0;h%solve_ITMAX(1)=100;h%solve_converged(1)=1
    endif
    if(all(rs==0)) then
     h%solve_SER(2)=0;h%solve_ITER(2)=0;h%solve_ITMAX(2)=100;h%solve_converged(2)=1
    endif
    h%worst_residual=maxval(h%solve_SER)
    h%cpcp(1)=1.e-3_real64*(maxval(h%pn)-minval(h%pn))
    h%cpcp(2)=1.e-3_real64*(maxval(h%ps)-minval(h%ps))
    h%facmax=1.e6_real64*max(maxval(abs(rn)),maxval(abs(rs)))
    call lfm_cube_yy_install_caps(h%pn,h%ps,rn,rs,step,time,h%spatial,local)
    write(9,*) 'YY_MI_STAGE INSTALL',step,time,local;flush(9)
   endif
  endif
  h%solve_seconds=h%solve_seconds+MPI_Wtime()-started
  call MPI_Bcast(local,1,MPI_INTEGER,0,h%comm,st);status=local;if(status/=0) return
  call MPI_Bcast(h%pn,size(h%pn),MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%ps,size(h%ps),MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%table,size(h%table),MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%spatial,size(h%spatial),MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%cpcp,2,MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%facmax,1,MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%worst_residual,1,MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%solve_SER,2,MPI_DOUBLE_PRECISION,0,h%comm,st)
  call MPI_Bcast(h%solve_ITER,2,MPI_INTEGER,0,h%comm,st)
  call MPI_Bcast(h%solve_ITMAX,2,MPI_INTEGER,0,h%comm,st)
  call MPI_Bcast(h%solve_converged,2,MPI_INTEGER,0,h%comm,st)
  call MPI_Bcast(h%mapped_trace_status,size(h%mapped_trace_status),MPI_INTEGER,0,h%comm,st)
  call yy_boundary_install(h%cap_x,h%cap_y,h%pn,h%ps,time,local)
  call collective_status(h%comm,local,status);if(status/=0) return
  started=MPI_Wtime()
  call yy_boundary_collective_refresh(h%comm,local)
  call collective_status(h%comm,local,status);if(status/=0) return
  call evaluate_wall(h,local)
  call collective_status(h%comm,local,status);if(status/=0) return
  h%wall_seconds=h%wall_seconds+MPI_Wtime()-started
  h%last_step=step;h%last_time=time;h%native_primed=.true.
  if(present(cpcp)) cpcp=h%cpcp
  if(present(facmax)) facmax=h%facmax
 end subroutine
 subroutine ion_wall_values(h,phi,status,time_s)
  class(cube_ionosphere),intent(in)::h
  real(real64),allocatable,intent(out)::phi(:,:)
  integer,intent(out)::status
  integer::j,k,st,node
  real(real64),optional,intent(in)::time_s
  status=1;if(.not.h%ready.or.h%last_step<0) return
  if(present(time_s)) then
   call yy_boundary_set_time(time_s,status);if(status/=0)return
  endif
  allocate(phi(0:h%layout%local(2),0:h%layout%local(3)));phi=0
  if(h%layout%coord(1)==0) then
   do k=0,h%layout%local(3)
    do j=0,h%layout%local(2)
     node=h%vertex(j+h%layout%first(2),k+h%layout%first(3),h%layout%face)
     if(present(time_s)) then
      ! Evaluate the same held magneticPhi at CURRENT edge endpoint geometry.
      ! Shared vertices use the same unique coordinate, no new solve/collective.
      call yy_boundary_potential(.01_real64*h%wall_xyz(:,node),phi(j,k),st)
      if(st/=0)return
     else
      phi(j,k)=h%wall_phi(node)
     endif
    enddo
   enddo
  endif
  status=0
 end subroutine
 subroutine ion_snapshot(h,data,status)
  class(cube_ionosphere),intent(in)::h
  real(real64),allocatable,intent(out)::data(:,:,:,:)
  integer,intent(out)::status
  status=1;if(.not.h%ready.or.h%last_step<0) return
  data=h%spatial;status=0
 end subroutine
 subroutine ion_export(h,buffer,status)
  class(cube_ionosphere),intent(in)::h
  real(real64),allocatable,intent(out)::buffer(:)
  integer,intent(out)::status
  call h%checkpoint_export(buffer,status)
 end subroutine
 subroutine ion_import(h,buffer,status)
  class(cube_ionosphere),intent(inout)::h
  real(real64),intent(in)::buffer(:)
  integer,intent(out)::status
  call h%checkpoint_import(buffer,status)
 end subroutine
 subroutine ion_checkpoint_export(h,buffer,status)
  class(cube_ionosphere),intent(in)::h
  real(real64),allocatable,intent(out)::buffer(:)
  integer,intent(out)::status
  character(64),parameter::coeff='e07b15a704841612ae7f0787e365bc75da46cb969720fd9e53143b9df9bd9bcd'
  integer::i
  status=1;if(.not.h%ready.or..not.h%native_primed.or.h%last_step<0) return
  buffer=[5._real64,100._real64,1.e-3_real64,1._real64,real(h%last_step,real64),h%last_time, &
   h%worst_residual,h%facmax,h%solve_SER,real(h%solve_ITER,real64),real(h%solve_converged,real64), &
   6380000._real64,rion,yy_spin_omega_rad_s,yy_spin_phase0_rad, &
   (real(iachar(coeff(i:i)),real64),i=1,64), &
   reshape(h%pn,[capj*capk]),reshape(h%ps,[capj*capk]), &
   reshape(h%spatial(:,:,6,1),[capj*capk]),reshape(h%spatial(:,:,6,2),[capj*capk]), &
   real(reshape(h%mapped_trace_status,[2*capj*capk]),real64)]
  if(.not.all(ieee_is_finite(buffer))) return
  status=0
 end subroutine
 subroutine ion_checkpoint_import(h,buffer,status)
  class(cube_ionosphere),intent(inout)::h
  real(real64),intent(in)::buffer(:)
  integer,intent(out)::status
  real(real64)::cn(capj,capk),cs(capj,capk)
  character(64),parameter::coeff='e07b15a704841612ae7f0787e365bc75da46cb969720fd9e53143b9df9bd9bcd'
  integer::n,local,st,i
  status=1
  if(.not.h%ready.or..not.h%native_primed.or.size(buffer)/=82+6*capj*capk) return
  if(.not.all(ieee_is_finite(buffer))) return
  if(any(buffer(:4)/=[5._real64,100._real64,1.e-3_real64,1._real64])) return
  if(any(buffer(15:18)/=[6380000._real64,rion,yy_spin_omega_rad_s,yy_spin_phase0_rad])) return
  if(any(buffer(19:82)/=[(real(iachar(coeff(i:i)),real64),i=1,64)])) return
  if(buffer(5)<0.or.buffer(5)/=anint(buffer(5)).or.buffer(6)<0) return
  if(any(buffer(9:10)<0).or.any(buffer(11:12)<0)) return
  if(any(buffer(13:14)/=real(merge(1,0,buffer(9:10)<=1.e-3_real64),real64))) return
  h%last_step=nint(buffer(5));h%last_time=buffer(6);h%worst_residual=buffer(7);h%facmax=buffer(8)
  h%solve_SER=buffer(9:10);h%solve_ITER=nint(buffer(11:12));h%solve_converged=nint(buffer(13:14))
  h%solve_ITMAX=100;n=82
  h%pn=reshape(buffer(n+1:n+capj*capk),shape(h%pn));n=n+capj*capk
  h%ps=reshape(buffer(n+1:n+capj*capk),shape(h%ps));n=n+capj*capk
  cn=reshape(buffer(n+1:n+capj*capk),shape(cn));n=n+capj*capk
  cs=reshape(buffer(n+1:n+capj*capk),shape(cs));n=n+capj*capk
  h%mapped_trace_status=reshape(nint(buffer(n+1:)),shape(h%mapped_trace_status))
  local=0
  if(h%rank==0) call lfm_cube_yy_install_caps(h%pn,h%ps,cn,cs,h%last_step,h%last_time,h%spatial,local)
  call MPI_Bcast(local,1,MPI_INTEGER,0,h%comm,st);status=local;if(status/=0) return
  call MPI_Bcast(h%spatial,size(h%spatial),MPI_DOUBLE_PRECISION,0,h%comm,st)
  h%cpcp(1)=1.e-3_real64*(maxval(h%pn)-minval(h%pn))
  h%cpcp(2)=1.e-3_real64*(maxval(h%ps)-minval(h%ps))
  call yy_boundary_install(h%cap_x,h%cap_y,h%pn,h%ps,h%last_time,local)
  call collective_status(h%comm,local,status);if(status/=0) return
  call yy_boundary_collective_refresh(h%comm,local)
  call collective_status(h%comm,local,status);if(status/=0) return
  call evaluate_wall(h,local);call collective_status(h%comm,local,status)
 end subroutine
end module
