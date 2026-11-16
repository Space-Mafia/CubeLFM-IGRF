! Reuse released trace/factor and interpolation on actual Cube source support.
module lfm_cube_yy_fac_map
 use iso_fortran_env,only:real64,error_unit
 use iso_c_binding,only:c_int,c_double
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 use lfm_cube_yy_native_cap_frame,only:yy_native_cap_frame
 use lfm_cube_mesh
 use lfm_cube_topology
 use lfm_cube_yy_mapping,only:yy_shell_intersection
 use lfm_cube_yy_igrf,only:cube_yy_igrf_field
 implicit none
 private
 public::yy_released_source_map
 interface
  integer(c_int)function interpolate(logical_xyz,face,n,p,value,clipped) &
    bind(C,name='cube_released_ghost_scalar_interpolate')
   import c_int,c_double
   real(c_double),intent(in)::logical_xyz(3),p(*)
   integer(c_int),value::face,n
   real(c_double),intent(out)::value
   integer(c_int),intent(out)::clipped
  end function
 end interface
contains
 subroutine yy_released_source_map(mesh,length_m,global_parallel,x,y,theta,radii,rion_m,time, &
  rn,rs,trace,status,work_first,work_last)
  type(cube_meridian),intent(in)::mesh
  real(real64),intent(in)::length_m
  real(real64),contiguous,intent(in)::global_parallel(0:,0:,0:,0:)
  real(real64),intent(in)::x(:,:),y(:,:),theta(:),radii(:),rion_m,time
  real(real64),intent(out)::rn(:,:),rs(:,:)
  integer,intent(out)::trace(:,:,:),status
  integer,intent(in)::work_first,work_last
  integer::j,k,h,flat,st,f,clip,nang
  real(real64)::t,p,foot(3),point(3),reference(3),a,b,logical_xyz(3),z,bmag,imag,direction(3),inc,q(3,3)
  real(real64),parameter::pi=acos(-1._real64)
  rn=0;rs=0;trace=-999;status=1
  q=yy_native_cap_frame(time)
  nang=size(global_parallel,2)-2
  if(length_m<=0.or.nang<4.or.size(global_parallel,1)/=6.or.size(global_parallel,4)/=6)return
  do h=1,2
   do k=3,size(x,2)-2
    do j=2,size(x,1)-1
     flat=(h-1)*size(x)+(k-1)*size(x,1)+j
     if(flat<work_first.or.flat>work_last)cycle
     t=hypot(x(j,k),y(j,k));p=modulo(atan2(y(j,k),x(j,k)),2*pi)
     if(h==2)t=pi-t
     foot=rion_m*matmul(q,[sin(t)*cos(p),sin(t)*sin(p),cos(t)])
     call yy_shell_intersection(foot,time,theta,radii,point,st)
     trace(j,k,h)=st
     ! Exact released policy: failed upward trace zero, status retained.
     if(st/=0)cycle
     call mesh%reference(point*100._real64,t,reference,st)
     if(st/=0)then
      write(error_unit,*)'YY_FAC_SOURCE_REFERENCE_FAILURE',h,j,k,st,point
      return
     endif
     call cube_inverse(reference,f,a,b,st)
     if(st/=0)return
     logical_xyz=[t*mesh%ni+.5_real64,(a/(pi/4)+1)*nang/2+.5_real64, &
      (b/(pi/4)+1)*nang/2+.5_real64]
     st=interpolate(logical_xyz,int(f,c_int),int(nang,c_int),global_parallel,z,clip)
     if(st/=0.or.clip/=0)then
      write(error_unit,*)'YY_FAC_SOURCE_STENCIL_FAILURE',h,j,k,st,clip,f,logical_xyz,point
      return
     endif
     st=cube_yy_igrf_field(point/length_m,time,direction,bmag)
     if(st/=0)return
     ! Same released signed inclination evaluated at the ionospheric foot.
     t=hypot(x(j,k),y(j,k));p=modulo(atan2(y(j,k),x(j,k)),2*pi)
     if(h==2)t=pi-t
     foot=matmul(q,[sin(t)*cos(p),sin(t)*sin(p),cos(t)])
     st=cube_yy_igrf_field(rion_m/length_m*foot,time,direction,imag)
     if(st/=0)return
     inc=dot_product(direction,foot);z=z*(imag/bmag)*inc
     if(h==1)rn(j,k)=z
     if(h==2)rs(j,k)=z
    enddo
   enddo
  enddo
  if(.not.all(ieee_is_finite(rn)).or..not.all(ieee_is_finite(rs)))return
  status=0
 end subroutine
end module
