module cube_event_global
 use iso_c_binding
 use iso_fortran_env, only: real64
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 use mpi_f08
 use cube_event_aurora_cap, only: event_np,event_channels,auroral_point,feet,advance_dpb,masks,precip
 use lfm_cube_yy_native_cap_frame, only: yy_native_cap_frame
 use lfm_cube_yy_fac_source, only: yy_source_global_moments
 implicit none
 private
 public :: global_update,global_caps,global_fields,global_metrics,global_held,global_size
 public :: nlon,nlat,global_write
 ! Two MI samples per angular MHD cell around four 24-cell faces.
 integer,parameter :: nlon=192,nlat=nlon/2+1,ng=nlon*nlat,global_size=4*ng
 real(real64),parameter :: pi=acos(-1._real64)
 include 'yy_rotation_profile.inc'
 ! GEO body-fixed east longitude, north-to-south colatitude, periodic without a duplicate seam.
 ! Channels: potential (V), P/H (S), outward radial FAC (A/m^2).
 real(real64),save :: global_fields(nlon,nlat,4)=0,global_metrics(4)=0,global_held=-1
 real(real64),save :: diagnostics(nlon,nlat,event_channels)=0
 interface
  integer(c_int) function solve(np,nt,r,f,p,h,v,m,it) bind(C,name='cube_event_global_solve')
   import
   integer(c_int),value :: np,nt
   real(c_double),value :: r
   real(c_double),intent(in) :: f(*),p(*),h(*)
   real(c_double),intent(inout) :: v(*)
   real(c_double),intent(out) :: m(*)
   integer(c_int),intent(out) :: it
  end function
 end interface
contains
 pure real(real64) function sample_grid(a,theta,east) result(v)
  real(real64),intent(in)::a(nlon,nlat),theta,east
  real(real64)::p,t,fp,ft
  integer::j,k,i
  p=modulo(east,2*pi)*(nlon/(2*pi))
  j=modulo(floor(p),nlon)+1;k=modulo(j,nlon)+1;fp=p-floor(p)
  t=min(real(nlat-1,real64),max(0._real64,theta*(real(nlat-1,real64)/pi)))
  i=min(nlat-1,int(t)+1);ft=t-(i-1)
  v=(1-ft)*((1-fp)*a(j,i)+fp*a(k,i))+ft*((1-fp)*a(j,i+1)+fp*a(k,i+1))
 end function

 integer(c_int) function return_potential(th,ph,value) bind(C,name='__wrap_cube_yy_footpoint_potential') result(st)
  real(c_double),value :: th,ph
  real(c_double),intent(out)::value
  st=1;value=0
  if(global_held<0.or..not.ieee_is_finite(th).or..not.ieee_is_finite(ph))return
  value=sample_grid(global_fields(:,:,1),th,ph-yy_spin_phase0_rad-yy_spin_omega_rad_s*global_held)
  if(ieee_is_finite(value))st=0
 end function

 subroutine global_update(comm,t,theta,radii,rion,ready,radius,edge,iterations,status)
  type(MPI_Comm),intent(in)::comm
  real(real64),intent(in)::t,theta(:),radii(:),rion
  integer,intent(inout)::ready(2)
  real(real64),intent(inout)::radius(2),edge(event_np,2)
  integer,intent(out)::iterations,status
  integer::rank,nproc,st,first,last,i,j,k,h,n,rc
  integer,allocatable :: counts(:),offsets(:)
  real(real64),allocatable::v(:,:),local(:,:),mom(:,:),z(:),cut(:),lat(:),mlt(:)
  real(real64),allocatable::dp(:,:,:,:),f(:,:,:),p(:,:,:),mask(:,:,:)
  integer,parameter::nmag=99
  integer::oldready(2)
  real(real64)::angle,th,ph,radial(nlon,nlat),prev_p,prev_h
  type(auroral_point)::point
  call MPI_Comm_rank(comm,rank,st);call MPI_Comm_size(comm,nproc,st)
  allocate(counts(0:nproc-1),offsets(0:nproc-1),v(3,ng));angle=yy_spin_phase0_rad+yy_spin_omega_rad_s*t
  do i=1,nlat
   th=(i-1)*pi/(nlat-1)
   do j=1,nlon
    k=j+(i-1)*nlon;ph=angle+(j-1)*(2*pi/nlon)
    v(:,k)=[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]
    if(i==1.or.i==nlat)v(:,k)=[0._real64,0._real64,sign(1._real64,cos(th))]
   enddo
  enddo
  do i=0,nproc-1
   offsets(i)=ng*i/nproc;counts(i)=ng*(i+1)/nproc-offsets(i)
  enddo
  first=offsets(rank)+1;last=first+counts(rank)-1
  allocate(local(5,counts(rank)),mom(5,ng))
  call yy_source_global_moments(t,theta,radii,rion,v(:,first:last),local,rc)
  call MPI_Allreduce(rc,status,1,MPI_INTEGER,MPI_MAX,comm,st)
  if(st/=MPI_SUCCESS.or.status/=0)then;status=1;return;endif
  call MPI_Gatherv(local,size(local),MPI_DOUBLE_PRECISION,mom,5*counts(:nproc-1),5*offsets(:nproc-1),MPI_DOUBLE_PRECISION,0,comm,st)
  if(st/=MPI_SUCCESS)then;status=1;return;endif
  iterations=0
  if(rank==0)then
   radial=reshape(mom(3,:)*mom(4,:),shape(radial))
   allocate(dp(3,event_np,nmag,2),f(event_np,nmag,2),p(event_np,nmag,2),mask(event_np,nmag,2))
   call feet(nmag,event_np,t,dp)
   do h=1,2
    do i=1,nmag
     do j=1,event_np
      th=atan2(norm2(dp(1:2,j,i,h)),dp(3,j,i,h));ph=atan2(dp(2,j,i,h),dp(1,j,i,h))-angle
      f(j,i,h)=sample_grid(radial,th,ph)
      p(j,i,h)=sample_grid(global_fields(:,:,1),th,ph)
     enddo
    enddo
   enddo
   oldready=ready
   rc=advance_dpb(nmag,event_np,f,p,ready,radius,edge,mask)
   allocate(cut(ng),z(ng),lat(ng),mlt(ng))
   if(rc==0)rc=masks(ng,event_np,t,v,edge,cut,z,lat,mlt)
   do k=1,ng
    if(rc/=0)exit
    i=1+(k-1)/nlon;j=1+modulo(k-1,nlon);h=merge(1,2,lat(k)>=0)
    prev_p=global_fields(j,i,2);prev_h=global_fields(j,i,3)
    rc=precip(t,z(k),merge(1,0,mom(5,k)==0),mom(1,k),mom(2,k),mom(3,k),mom(4,k),cut(k),oldready(h),prev_p,prev_h,point)
    if(rc/=0)exit
    global_fields(j,i,2:3)=[point%sp,point%sh]
    diagnostics(j,i,:)=[point%energy,point%number,point%mean,point%ep,point%eh,point%ap,point%ah,cut(k),z(k),real(point%diffuse,real64),lat(k),mlt(k)]
   enddo
   global_fields(:,:,4)=radial
   if(rc==0)rc=solve(nlon,nlat,rion,radial,global_fields(:,:,2),global_fields(:,:,3),global_fields(:,:,1),global_metrics,iterations)
  endif
  call MPI_Bcast(rc,1,MPI_INTEGER,0,comm,st)
  if(st/=MPI_SUCCESS.or.rc/=0)then;status=1;return;endif
  call MPI_Bcast(global_fields,size(global_fields),MPI_DOUBLE_PRECISION,0,comm,st)
  if(st/=MPI_SUCCESS)then;status=1;return;endif
  call MPI_Bcast(global_metrics,4,MPI_DOUBLE_PRECISION,0,comm,st)
  if(st/=MPI_SUCCESS)then;status=1;return;endif
  global_held=t;status=0
 end subroutine

 subroutine global_caps(t,x,y,pn,ps,sp,sh,diag)
  real(real64),intent(in)::t,x(:,:),y(:,:)
  real(real64),intent(out)::pn(:,:),ps(:,:),sp(:,:,:),sh(:,:,:),diag(:,:,:,:)
  real(real64)::q(3,3),v(3),z,a,b,th,lon,phi
  integer::i,j,h,k
  q=yy_native_cap_frame(t)
  do h=1,2
   do j=1,size(x,2)
    do i=1,size(x,1)
     a=hypot(x(i,j),y(i,j));b=atan2(y(i,j),x(i,j));z=cos(a)
     if(h==2)z=-z
     v=matmul(q,[sin(a)*cos(b),sin(a)*sin(b),z])
     th=atan2(norm2(v(:2)),v(3));lon=atan2(v(2),v(1))-yy_spin_phase0_rad-yy_spin_omega_rad_s*t
     phi=sample_grid(global_fields(:,:,1),th,lon)
     if(h==1)pn(i,j)=phi
     if(h==2)ps(i,j)=phi
     sp(i,j,h)=sample_grid(global_fields(:,:,2),th,lon);sh(i,j,h)=sample_grid(global_fields(:,:,3),th,lon)
     do k=1,event_channels
      diag(k,i,j,h)=sample_grid(diagnostics(:,:,k),th,lon)
     enddo
     ! Preserve categorical diffuse flags and avoid interpolating across the MLT wrap.
     block
      integer::jj,ii
      jj=1+modulo(nint(modulo(lon,2*pi)*nlon/(2*pi)),nlon)
      ii=1+min(nlat-1,max(0,nint(th*(nlat-1)/pi)))
      diag(10,i,j,h)=diagnostics(jj,ii,10);diag(12,i,j,h)=diagnostics(jj,ii,12)
     end block
    enddo
   enddo
  enddo
 end subroutine

 subroutine global_write(prefix,step,time)
  character(*),intent(in)::prefix
  integer,intent(in)::step
  real(real64),intent(in)::time
  integer::u
  open(newunit=u,file=trim(prefix)//'.io-tmp/global-state',status='new',access='stream',form='unformatted')
  write(u) step,nlon,nlat,time,global_held,global_metrics,global_fields
  close(u)
 end subroutine
end module
