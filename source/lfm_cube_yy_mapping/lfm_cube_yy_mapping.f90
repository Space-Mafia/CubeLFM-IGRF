! Layout/units adapter: qualified YY tracer and original LFM interpolation.
! No new field-line integrator, elliptic operator or extrapolation scheme.
module lfm_cube_yy_mapping
 use iso_fortran_env,only:real64
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 use lfm_cube_yy_igrf,only:cube_yy_igrf_configure,cube_yy_igrf_field,cube_yy_igrf_trace
 use lfm_rotating_mi_map,only:mi_cap_weights,mi_shell_scalar
 implicit none
 private
 public::yy_mapping_configure,yy_map_fac,yy_cap_potential,yy_corotation_edge,yy_shell_intersection
 real(real64),parameter::pi=acos(-1._real64),field_unit_t=1.e-9_real64
 real(real64),save::length_m=0,rion_m=0,omega=0,phase0=0
 logical,save::ready=.false.
 ! Native GINT callback context; this branch requires OMP1 per MPI process.
 real(real64),save::corotation_epoch=0,corotation_edge_vector(3)=0
 logical,save::corotation_failed=.false.
contains
 subroutine yy_mapping_configure(model_length_m,ionosphere_m,omega_s,phase_rad,status)
  real(real64),intent(in)::model_length_m,ionosphere_m,omega_s,phase_rad
  integer,intent(out)::status
  status=1
  if(.not.all(ieee_is_finite([model_length_m,ionosphere_m,omega_s,phase_rad]))) return
  if(min(model_length_m,ionosphere_m)<=0) return
  status=cube_yy_igrf_configure(model_length_m,field_unit_t,omega_s,phase_rad)
  if(status/=0) return
  length_m=model_length_m;rion_m=ionosphere_m;omega=omega_s;phase0=phase_rad;ready=.true.
 end subroutine

 subroutine yy_map_fac(x,y,theta,phi,radii,jcart,time_s,smooth,rn,rs,trace_status,status,work_first,work_last)
  ! Native JPARA Cartesian curl is supplied in A/m^2. Its scalar projection
  ! and the signed ionospheric B ratio follow released fill_local_source and
  ! map_mhd_source_to_ionosphere. Native angular interpolation is reused.
  real(real64),intent(in)::x(:,:),y(:,:),theta(:),phi(:),radii(:),jcart(:,:,:),time_s
  logical,intent(in)::smooth
  real(real64),intent(out)::rn(:,:),rs(:,:)
  integer,intent(out)::trace_status(:,:,:),status
  real(real64)::parallel(size(theta),size(phi)),filtered(size(theta),size(phi))
  real(real64)::point(3),foot(3),start(3),endpoint(4),direction(3),imag,bmag,z,inc,t,p,rsample
  integer,optional,intent(in)::work_first,work_last
  integer::j,k,h,st,n,m,flat,first,last
  status=1;rn=0;rs=0;trace_status=-999
  n=size(theta);m=size(phi)
  first=1;last=2*size(x)
  if(present(work_first)) first=work_first
  if(present(work_last)) last=work_last
  if(first<1.or.last>2*size(x).or.last<first-1) return
  if(.not.ready.or.n<3.or.m<5.or.size(radii)/=n) return
  if(any(shape(jcart)/=[3,n,m]).or.any(shape(x)/=shape(y))) return
  if(any(shape(rn)/=shape(x)).or.any(shape(rs)/=shape(x))) return
  if(any(shape(trace_status)/=[size(x,1),size(x,2),2])) return
  if(.not.all(ieee_is_finite(jcart)).or..not.all(ieee_is_finite(radii))) return
  ! Preserve the actual native shell. Only the tracing point provider changes;
  ! its existing radius(theta) intersection algorithm is reused below.
  if(minval(radii)<=rion_m) return
  do k=1,m
   do j=1,n
    point=radii(j)/length_m*[cos(theta(j)),sin(theta(j))*cos(phi(k)),sin(theta(j))*sin(phi(k))]
    st=cube_yy_igrf_field(point,time_s,direction,bmag)
    if(st/=0) return
    parallel(j,k)=dot_product(jcart(:,j,k),direction)
   enddo
  enddo
  if(smooth) then
   ! Exact existing LFM nine-point smoothing body, including its extents.
   filtered=parallel
   do k=2,m-2
    do j=2,n-1
     filtered(j,k)=.25_real64*parallel(j,k)+.125_real64*(parallel(j+1,k)+parallel(j,k+1) &
      +parallel(j-1,k)+parallel(j,k-1))+.0625_real64*(parallel(j-1,k-1)+parallel(j-1,k+1) &
      +parallel(j+1,k-1)+parallel(j+1,k+1))
    enddo
   enddo
   parallel=filtered
  endif
  do h=1,2
   do k=3,size(x,2)-2
    do j=2,size(x,1)-1
     flat=(h-1)*size(x)+(k-1)*size(x,1)+j
     if(flat<first.or.flat>last) cycle
     t=hypot(x(j,k),y(j,k));p=modulo(atan2(y(j,k),x(j,k)),2*pi)
     if(h==2) t=pi-t
     foot=rion_m*[sin(t)*cos(p),sin(t)*sin(p),cos(t)]
     call yy_shell_intersection(foot,time_s,theta,radii,point,st)
     trace_status(j,k,h)=st
     ! Released upward trace failure: zero FAC, with the failure retained.
     ! It is not classified as a successful zero-current mapping.
     if(st/=0) cycle
     point=point/length_m
     call mi_shell_scalar(point,theta,phi,parallel,z,st)
     if(st/=0) then
      status=3;return
     endif
     st=cube_yy_igrf_field(point,time_s,direction,bmag)
     if(st/=0) return
     foot=[sin(t)*cos(p),sin(t)*sin(p),cos(t)]
     st=cube_yy_igrf_field(rion_m/length_m*foot,time_s,direction,imag)
     if(st/=0) return
     inc=dot_product(direction,foot)
     z=z*(imag/bmag)*inc
     if(h==1) rn(j,k)=z
     if(h==2) rs(j,k)=z
    enddo
   enddo
  enddo
  if(.not.all(ieee_is_finite(rn)).or..not.all(ieee_is_finite(rs))) return
  status=0
 end subroutine

 subroutine yy_cap_potential(point_m,time_s,x,y,pn,ps,value,trace_status,status)
  ! Qualified YY single footpoint plus unchanged native cap interpolation.
  ! No analytic dipole footpoint/Jacobian is used for IGRF.
  real(real64),intent(in)::point_m(3),time_s,x(:,:),y(:,:),pn(:,:),ps(:,:)
  real(real64),intent(out)::value
  integer,intent(out)::trace_status,status
  real(real64)::r,start(3),endpoint(4),t,p,w(4),dw(2,4),v(4)
  integer::ij(2),j,k
  value=0;status=1;trace_status=-999
  if(.not.ready.or..not.all(ieee_is_finite(point_m))) return
  if(any(shape(pn)/=shape(x)).or.any(shape(ps)/=shape(x))) return
  r=norm2(point_m);if(r<=0) return
  start=[r/length_m,acos(max(-1._real64,min(1._real64,point_m(3)/r))), &
   modulo(atan2(point_m(2),point_m(1)),2*pi)]
  trace_status=cube_yy_igrf_trace(start,rion_m/length_m,time_s,endpoint)
  if(trace_status/=0) return
  t=min(endpoint(1),pi-endpoint(1));p=endpoint(2)
  call mi_cap_weights(x,y,t*[cos(p),sin(p)],ij,w,dw,status)
  if(status/=0) return
  j=ij(1);k=ij(2)
  if(endpoint(1)<pi/2) then
   v=[pn(j,k),pn(j+1,k),pn(j,k+1),pn(j+1,k+1)]
  else
   v=[ps(j,k),ps(j+1,k),ps(j,k+1),ps(j+1,k+1)]
  endif
  value=dot_product(w,v)
  if(.not.ieee_is_finite(value)) status=1
 end subroutine

 subroutine yy_corotation_edge(a_m,b_m,time_s,emf_lfm,status)
  ! Same released v_rot x fullIGRF. Use the existing original LFM GINT
  ! instead of a midpoint sample on a general mapped-Cube straight edge.
  real(real64),intent(in)::a_m(3),b_m(3),time_s
  real(real64),intent(out)::emf_lfm
  integer,intent(out)::status
  real(real64),external::gint
  emf_lfm=0;status=1
  if(.not.ready.or..not.all(ieee_is_finite([a_m,b_m,time_s])))return
  corotation_epoch=time_s;corotation_edge_vector=b_m-a_m;corotation_failed=.false.
  emf_lfm=gint(yy_corotation_integrand,a_m(1),a_m(2),a_m(3),b_m(1),b_m(2),b_m(3))
  if(corotation_failed.or..not.ieee_is_finite(emf_lfm))return
  status=0
 end subroutine
 real(real64) function yy_corotation_integrand(x,y,z) result(value)
  real(real64),intent(in)::x,y,z
  real(real64)::point(3),direction(3),mag,b(3),v(3),vxB(3)
  integer::st
  value=0;if(corotation_failed)return
  point=[x,y,z]
  st=cube_yy_igrf_field(point/length_m,corotation_epoch,direction,mag)
  if(st/=0)then
   corotation_failed=.true.;return
  endif
  b=direction*mag*field_unit_t;v=[-omega*point(2),omega*point(1),0._real64]
  vxB=[v(2)*b(3)-v(3)*b(2),v(3)*b(1)-v(1)*b(3),v(1)*b(2)-v(2)*b(1)]
  value=1.e8_real64*dot_product(vxB,corotation_edge_vector)
  if(.not.ieee_is_finite(value))corotation_failed=.true.
 end function
  subroutine yy_shell_intersection(foot,time_s,theta,radii,r,status)
    ! Reuse the original native shell interpolation and bisection body.
    ! Replace only its analytic dipole point provider with qualified YY RK4.
    ! The native shell need not be a sphere. Coordinates and radii are SI.
    real(real64),intent(in)::foot(3),time_s,theta(:),radii(:)
    real(real64),intent(out)::r(3)
    integer,intent(out)::status
    real(real64)::lo,hi,mid,flo,fhi,fmid,tmp(3)
    integer::n,iter,st
    logical::bracket_found
    r=0;status=1;n=size(theta)
    if(n<2.or.size(radii)/=n.or.rion_m<=0) return
    if(.not.all(ieee_is_finite(theta)).or..not.all(ieee_is_finite(radii))) return
    if(any(theta(2:)<=theta(:n-1)).or.minval(radii)<=rion_m) return
    if(theta(1)<0.or.theta(1)>1.e-8_real64.or.theta(n)<pi-1.e-8_real64.or.theta(n)>pi) return
    if(abs(norm2(foot)/rion_m-1)>1.e-10_real64) return
    lo=minval(radii)*(1-1.e-12_real64);hi=maxval(radii)*(1+1.e-12_real64)
    if(hi<lo) then
      status=3;return
    endif
    call residual(lo,flo,tmp,st)
    if(st/=0) return
    call residual(hi,fhi,tmp,st)
    if(st/=0) then
      ! The line can intersect this angle-dependent LFM shell without
      ! reaching its largest radius anywhere on the shell. A failed trace
      ! to that bounding sphere therefore does not prove missing FAC.
      ! Reuse the same qualified YY tracer and native bisection: find a
      ! reachable nonnegative residual before solving the shell crossing.
      bracket_found=.false.
      do iter=1,60
        mid=(lo+hi)/2
        call residual(mid,fmid,tmp,st)
        if(st/=0) then
          hi=mid
        else
          if(abs(fmid)<2.e-14_real64*mid) then
            r=tmp;status=0;return
          endif
          if(fmid>=0) then
            hi=mid;fhi=fmid;bracket_found=.true.;exit
          endif
          lo=mid;flo=fmid
        endif
      enddo
      if(.not.bracket_found) return
    endif
    if(flo*fhi>0) then
      status=3;return
    endif
    do iter=1,60
      mid=(lo+hi)/2
      call residual(mid,fmid,tmp,st)
      if(st/=0) return
      if(abs(fmid)<2.e-14_real64*mid.or.hi-lo<4.e-14_real64*mid) then
        r=tmp;status=0;return
      endif
      if(flo*fmid<=0) then
        hi=mid
      else
        lo=mid;flo=fmid
      endif
    enddo
    status=4
  contains
    subroutine residual(rad,value,pos,ierr)
      real(real64),intent(in)::rad
      real(real64),intent(out)::value,pos(3)
      integer,intent(out)::ierr
      real(real64)::angle,u,start(3),endpoint(4)
      integer::j
      start=[rion_m/length_m,acos(max(-1._real64,min(1._real64,foot(3)/rion_m))), &
        modulo(atan2(foot(2),foot(1)),2*pi)]
      ierr=cube_yy_igrf_trace(start,rad/length_m,time_s,endpoint)
      value=0
      if(ierr/=0) return
      pos=rad*[sin(endpoint(1))*cos(endpoint(2)),sin(endpoint(1))*sin(endpoint(2)),cos(endpoint(1))]
      angle=atan2(norm2(pos(2:)),pos(1))
      j=1
      do while(j<n-1)
        if(theta(j+1)>=angle) exit
        j=j+1
      enddo
      u=max(0._real64,min(1._real64,(angle-theta(j))/(theta(j+1)-theta(j))))
      value=rad-((1-u)*radii(j)+u*radii(j+1))
    end subroutine
  end subroutine

end module
