! Fortran port of GAMERA-Cube8531526's LFM logical-shell geometry only.
! Coefficients, PCHIP, extension and inverse follow cubedsphere_mesh.c;
! target convention is LFM +X sunward (proper rotation of the C convention).
module lfm_cube_shell_map
 use iso_fortran_env,only:real64
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::lfm_shell_map
 type lfm_shell_map
  integer::ns=0,nt=0
  logical::direct_shells=.false.
  real(real64)::smin=0,s6=0,s10=0,rin=0,sphere=0,transition=0,extra=0
  real(real64),allocatable::s(:),theta(:),coefficient(:,:,:,:),seed(:,:)
 contains
  procedure::load=>shell_load
  procedure::forward=>shell_forward
  procedure::inverse=>shell_inverse
  procedure::meridian=>shell_meridian
 end type
contains
 integer function interval(knots,x) result(i)
  real(real64),intent(in)::knots(0:),x
  integer::lo,hi,mid
  lo=0;hi=size(knots)-1
  do while(hi-lo>1)
   mid=(lo+hi)/2
   if(x<knots(mid)) then
    hi=mid
   else
    lo=mid
   endif
  enddo
  i=lo
 end function
 real(real64) function angular(m,row,theta,component) result(y)
  class(lfm_shell_map),intent(in)::m
  integer,intent(in)::row,component
  real(real64),intent(in)::theta
  real(real64)::t,x,c(4)
  integer::j
  t=modulo(theta-m%theta(0),m%theta(m%nt-1)-m%theta(0))+m%theta(0)
  j=interval(m%theta,t);c=m%coefficient(:,component,j,row);x=t-m%theta(j)
  y=((c(1)*x+c(2))*x+c(3))*x+c(4)
 end function
 real(real64) function slope(m,i,theta,component) result(v)
  class(lfm_shell_map),intent(in)::m
  integer,intent(in)::i,component
  real(real64),intent(in)::theta
  integer::a,b,c
  real(real64)::h0,h1,d0,d1,w0,w1,y
  if(i==0.or.i==m%ns-1) then
   if(i==0) then
    a=0;b=1;c=2
   else
    a=m%ns-1;b=m%ns-2;c=m%ns-3
   endif
   h0=abs(m%s(b)-m%s(a));h1=abs(m%s(c)-m%s(b))
   d0=(angular(m,b,theta,component)-angular(m,a,theta,component))/(m%s(b)-m%s(a))
   d1=(angular(m,c,theta,component)-angular(m,b,theta,component))/(m%s(c)-m%s(b))
   v=((2*h0+h1)*d0-h0*d1)/(h0+h1)
   if(v*d0<=0) then
    v=0
   else if(d0*d1<=0.and.abs(v)>3*abs(d0)) then
    v=3*d0
   endif
  else
   h0=m%s(i)-m%s(i-1);h1=m%s(i+1)-m%s(i);y=angular(m,i,theta,component)
   d0=(y-angular(m,i-1,theta,component))/h0;d1=(angular(m,i+1,theta,component)-y)/h1
   v=0
   if(d0*d1<=0) return
   w0=2*h1+h0;w1=h1+2*h0;v=(w0+w1)/(w0/d0+w1/d1)
  endif
 end function
 real(real64) function radial(m,s,theta,component,derivative) result(v)
  class(lfm_shell_map),intent(in)::m
  real(real64),intent(in)::s,theta
  integer,intent(in)::component
  real(real64),intent(out),optional::derivative
  integer::i
  real(real64)::h,x,y,z,a,b,c,d
  i=interval(m%s,s);h=m%s(i+1)-m%s(i);x=s-m%s(i)
  y=angular(m,i,theta,component);z=angular(m,i+1,theta,component)
  a=slope(m,i,theta,component);b=slope(m,i+1,theta,component)
  d=(z-y)/h;c=(3*d-2*a-b)/h;b=(a+b-2*d)/(h*h)
  if(present(derivative)) derivative=(3*b*x+2*c)*x+a
  v=((b*x+c)*x+a)*x+y
 end function
 real(real64) function smooth(x) result(v)
  real(real64),intent(in)::x
  real(real64)::y
  y=max(0._real64,min(1._real64,x));v=y*y*y*(10+y*(-15+6*y))
 end function
 real(real64) function radial_fraction(m,s) result(v)
  class(lfm_shell_map),intent(in)::m
  real(real64),intent(in)::s
  if(m%direct_shells) then
   v=s;return
  endif
  v=((s-m%smin)/(1-m%smin)+m%extra*smooth((s-m%s6)/(m%s10-m%s6)))/(1+m%extra)
 end function
 real(real64) function source_fraction(m,t) result(s)
  class(lfm_shell_map),intent(in)::m
  real(real64),intent(in)::t
  real(real64)::lo,hi,mid
  integer::j
  if(m%direct_shells) then
   s=t;return
  endif
  if(t<=0) then
   s=m%smin+t*(1-m%smin)*(1+m%extra);return
  else if(t>=1) then
   s=1+(t-1)*(1-m%smin)*(1+m%extra);return
  endif
  lo=m%smin;hi=1
  do j=1,digits(s)
   mid=(lo+hi)*.5_real64
   if(mid==lo.or.mid==hi) exit
   if(radial_fraction(m,mid)<t) then
    lo=mid
   else
    hi=mid
   endif
  enddo
  s=(lo+hi)*.5_real64
 end function
 subroutine shell_meridian(m,s,theta,x,r)
  class(lfm_shell_map),intent(in)::m
  real(real64),intent(in)::s,theta
  real(real64),intent(out)::x,r
  real(real64)::ref,derivative,w,dx,dr,factor
  if(m%direct_shells.and.s<=1) then
   if(s>=0) then
    x=radial(m,s,theta,0);r=radial(m,s,theta,1)
   else
    x=radial(m,0._real64,theta,0,dx);r=radial(m,0._real64,theta,1,dr)
    factor=exp(s*(x*dx+r*dr)/(x*x+r*r));x=x*factor;r=r*factor
   endif
   return
  endif
  if(s<m%smin) then
   ref=radial(m,m%smin,0._real64,0,derivative)
   ref=m%rin*exp((s-m%smin)*derivative/m%rin);x=ref*cos(theta);r=ref*sin(theta)
  else if(s>1) then
   x=radial(m,1._real64,theta,0,derivative);x=x+(s-1)*derivative
   r=radial(m,1._real64,theta,1,derivative);r=r+(s-1)*derivative
  else
   ref=radial(m,s,0._real64,0);w=smooth((ref-m%sphere)/(m%transition-m%sphere))
   x=ref*cos(theta)+w*(radial(m,s,theta,0)-ref*cos(theta))
   r=ref*sin(theta)+w*(radial(m,s,theta,1)-ref*sin(theta))
  endif
 end subroutine
 subroutine shell_load(m,path,status)
  class(lfm_shell_map),intent(out)::m
  character(*),intent(in)::path
  integer,intent(out)::status
  integer::u,st,i,j,k
  character(80)::magic
  real(real64)::values(8),x,r,t
  status=1
  open(newunit=u,file=path,status='old',action='read',iostat=st)
  if(st/=0) return
  read(u,*,iostat=st) magic
  if(st/=0.or.(magic/='GAMERA_CS_LFMD_V1'.and.magic/='LFM_CUBE_SHELL_V1')) then
   close(u);return
  endif
  m%direct_shells=magic=='LFM_CUBE_SHELL_V1'
  read(u,*,iostat=st) m%ns,m%nt
  if(st/=0.or.m%ns<3.or.m%nt<3.or.m%ns>10000.or.m%nt>10000) then
   close(u);return
  endif
  if(m%direct_shells) then
   read(u,*,iostat=st) m%rin
  else
   read(u,*,iostat=st) m%smin,m%s6,m%s10,m%rin,m%sphere,m%transition,m%extra
  endif
  if(st/=0) then
   close(u);return
  endif
  if(m%direct_shells) then
   if(m%rin<=0.or..not.ieee_is_finite(m%rin)) then
    close(u);return
   endif
  else if(.not.(m%smin>=0.and.m%smin<m%s6.and.m%s6<m%s10.and.m%s10<1.and. &
    m%rin>0.and.m%rin<m%sphere.and.m%sphere<m%transition.and.m%extra>=0)) then
   close(u);return
  endif
  allocate(m%s(0:m%ns-1),m%theta(0:m%nt-1),m%coefficient(4,0:1,0:m%nt-2,0:m%ns-1))
  read(u,*,iostat=st) m%s
  if(st==0) read(u,*,iostat=st) m%theta
  do i=0,m%ns-1
   do j=0,m%nt-2
    if(st/=0) exit
    read(u,*,iostat=st) values
    m%coefficient(:,0,j,i)=values(:4);m%coefficient(:,1,j,i)=values(5:)
   enddo
  enddo
  close(u)
  if(st/=0) return
  if(.not.all(ieee_is_finite(m%coefficient)).or..not.all(ieee_is_finite(m%s)).or. &
    .not.all(ieee_is_finite(m%theta))) return
  if(any(m%s(1:)<=m%s(:m%ns-2)).or.any(m%theta(1:)<=m%theta(:m%nt-2))) return
  if(m%s(0)/=0.or.m%s(m%ns-1)/=1.or.m%theta(0)/=0) return
  if(abs(m%theta(m%nt-1)-2*acos(-1._real64))>epsilon(t)*2*acos(-1._real64)) return
  allocate(m%seed(4,m%ns*((m%nt+1)/2)));k=0
  do i=0,m%ns-1
   do j=0,(m%nt-1)/2
    t=m%smin+(1-m%smin)*m%s(i)
    call m%meridian(t,m%theta(j),x,r);k=k+1;m%seed(:,k)=[t,m%theta(j),x,r]
   enddo
  enddo
  status=0
 end subroutine
 subroutine shell_forward(m,t,direction,point,status)
  class(lfm_shell_map),intent(in)::m
  real(real64),intent(in)::t,direction(3)
  real(real64),intent(out)::point(3)
  integer,intent(out)::status
  real(real64)::length,s,ref,transverse,theta,x,r
  status=1;point=0;length=norm2(direction)
  if(m%ns<3.or..not.ieee_is_finite(t).or.length<=0) return
  s=source_fraction(m,t);ref=0
  if(.not.m%direct_shells) ref=radial(m,s,0._real64,0)
  if(.not.m%direct_shells.and.s>=m%smin.and.ref<=m%sphere) then
   point=ref*direction/length
  else
   transverse=norm2(direction(2:));theta=atan2(transverse,direction(1))
   call m%meridian(s,theta,x,r);point(1)=x
   if(transverse>0) point(2:)=r*direction(2:)/transverse
  endif
  if(all(ieee_is_finite(point))) status=0
 end subroutine
 subroutine shell_inverse(m,point,t,direction,status)
  class(lfm_shell_map),intent(in)::m
  real(real64),intent(in)::point(3)
  real(real64),intent(out)::t,direction(3)
  integer,intent(out)::status
  real(real64)::x,r,length,s,theta,derivative,ref,lo,hi,mid,nearest,error,tolerance,h,hs,ht
  real(real64)::fx,fr,xp,rp,xm,rm,a,b,c,d,det,ds,dt,factor,ns,nt
  logical::converged,accepted
  integer::i,j,k
  status=1;t=0;direction=0;x=point(1);r=norm2(point(2:));length=norm2(point)
  if(m%ns<3.or.length<=0.or..not.all(ieee_is_finite(point))) return
  s=m%smin;theta=atan2(r,x)
  if(.not.m%direct_shells.and.length<=m%sphere) then
   if(length<m%rin) then
    ref=radial(m,m%smin,0._real64,0,derivative);s=m%smin+log(length/m%rin)*m%rin/derivative
   else
    lo=m%smin;hi=m%s6
    do j=1,digits(s)
     mid=(lo+hi)*.5_real64
     if(mid==lo.or.mid==hi) exit
     if(radial(m,mid,0._real64,0)<length) then
      lo=mid
     else
      hi=mid
     endif
    enddo
    s=(lo+hi)*.5_real64
   endif
   t=radial_fraction(m,s);direction=point/length;status=0;return
  endif
  nearest=huge(s)
  do i=1,size(m%seed,2)
   error=norm2([x-m%seed(3,i),r-m%seed(4,i)])
   if(error<nearest) then
    nearest=error;s=m%seed(1,i);theta=m%seed(2,i)
   endif
  enddo
  tolerance=256*epsilon(s)*max(1._real64,length);h=epsilon(s)**(1._real64/3);converged=.false.
  do i=1,digits(s)
   call m%meridian(s,theta,fx,fr);error=norm2([fx-x,fr-r])
   if(error<=tolerance) then
    converged=.true.;exit
   endif
   hs=h*max(1._real64,abs(s));ht=h*max(1._real64,abs(theta))
   call m%meridian(s+hs,theta,xp,rp);call m%meridian(s-hs,theta,xm,rm)
   a=(xp-xm)/(2*hs);c=(rp-rm)/(2*hs)
   call m%meridian(s,theta+ht,xp,rp);call m%meridian(s,theta-ht,xm,rm)
   b=(xp-xm)/(2*ht);d=(rp-rm)/(2*ht);det=a*d-b*c
   if(det==0.or..not.ieee_is_finite(det)) return
   ds=(d*(fx-x)-b*(fr-r))/det;dt=(a*(fr-r)-c*(fx-x))/det
   accepted=.false.;factor=1
   do k=1,digits(s)
    ns=s-factor*ds;nt=max(0._real64,min(acos(-1._real64),theta-factor*dt))
    call m%meridian(ns,nt,xp,rp)
    if(all(ieee_is_finite([xp,rp])).and.norm2([xp-x,rp-r])<error) then
     s=ns;theta=nt;accepted=.true.;exit
    endif
    factor=factor*.5_real64
   enddo
   if(.not.accepted) return
  enddo
  if(.not.converged) return
  t=radial_fraction(m,s);direction(1)=cos(theta)
  if(r>0) direction(2:)=sin(theta)*point(2:)/r
  status=0
 end subroutine
end module
