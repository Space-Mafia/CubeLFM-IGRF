! Stateless dynamic dipole MI geometry. SI units throughout this module.
! No MHD, FAC curl, or ionospheric elliptic discretization is replaced here.
module lfm_rotating_mi_geometry
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lfm_tilt_geometry, only: dipole_rotation, dipole_field, cross3
  implicit none
  private
  integer, parameter :: rk=real64
  real(rk), parameter :: pi=acos(-1._rk)
  public :: mi_frame, mi_footpoint_jacobian, mi_outward_point
  public :: mi_fac_radial, mi_pullback_e, mi_edge_lfm, mi_wall_drift
  public :: mi_shell_intersection
contains
  pure subroutine mi_shell_intersection(foot,q,rion,theta,radii,r,status)
    ! Intersect the outward dipole branch with an AXISYMMETRIC native shell.
    ! theta is colatitude about computational +X (not the magnetic axis).
    ! Shell model: linear interpolation of nodal radii versus theta, with
    ! axis values extended only over the grid's small endpoint epsilon.
    ! This solves one scalar equation, never numerically integrates a field
    ! line and never replaces a slightly non-spherical shell by its mean.
    ! status=3: no bracket before dipole apex (inaccessible/grazing footprint).
    real(rk),intent(in)::foot(3),q(3,3),rion,theta(:),radii(:)
    real(rk),intent(out)::r(3)
    integer,intent(out)::status
    real(rk)::fm(3),s2,apex,lo,hi,mid,flo,fhi,fmid,tmp(3)
    integer::n,iter,st
    r=0;status=1;n=size(theta)
    if(n<2.or.size(radii)/=n.or.rion<=0) return
    if(.not.all(ieee_is_finite(theta)).or..not.all(ieee_is_finite(radii))) return
    if(any(theta(2:)<=theta(:n-1)).or.minval(radii)<=rion) return
    if(theta(1)<0.or.theta(1)>1.e-8_rk.or.theta(n)<pi-1.e-8_rk.or.theta(n)>pi) return
    if(abs(norm2(foot)/rion-1)>1.e-10_rk) return
    fm=matmul(transpose(q),foot);s2=sum((fm(:2)/rion)**2)
    lo=minval(radii)*(1-1.e-12_rk);hi=maxval(radii)*(1+1.e-12_rk)
    if(s2>tiny(s2)) then
      apex=rion/s2
      hi=min(hi,apex)
    endif
    if(hi<lo) then
      status=3;return
    endif
    call residual(lo,flo,tmp,st)
    if(st/=0) return
    call residual(hi,fhi,tmp,st)
    if(st/=0) return
    if(flo*fhi>0) then
      status=3;return
    endif
    do iter=1,60
      mid=(lo+hi)/2
      call residual(mid,fmid,tmp,st)
      if(st/=0) return
      if(abs(fmid)<2.e-14_rk*mid.or.hi-lo<4.e-14_rk*mid) then
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
    pure subroutine residual(rad,value,pos,ierr)
      real(rk),intent(in)::rad
      real(rk),intent(out)::value,pos(3)
      integer,intent(out)::ierr
      real(rk)::angle,u
      integer::j
      call mi_outward_point(foot,q,rion,rad,pos,ierr)
      value=0
      if(ierr/=0) return
      angle=atan2(norm2(pos(2:)),pos(1))
      j=1
      do while(j<n-1)
        if(theta(j+1)>=angle) exit
        j=j+1
      enddo
      u=max(0._rk,min(1._rk,(angle-theta(j))/(theta(j+1)-theta(j))))
      value=rad-((1-u)*radii(j)+u*radii(j+1))
    end subroutine
  end subroutine

  pure function mi_frame(spin_frame,alpha,phase0,omega,time) result(q)
    real(rk), intent(in) :: spin_frame(3,3),alpha,phase0,omega,time
    real(rk) :: q(3,3)
    q=matmul(spin_frame,dipole_rotation(alpha,modulo(phase0+omega*time,2*pi)))
  end function

  pure subroutine mi_footpoint_jacobian(r,q,rion,foot,df,hemisphere,status)
    ! q: magnetic -> computational frame. df(a,b)=d foot(a)/d r(b).
    ! Each hemisphere is a separate inward branch. At magnetic equator
    ! the branch is ambiguous; status=2 requires a cap-rim policy, not a
    ! silently chosen magnetic hemisphere. Exact magnetic poles are regular.
    real(rk), intent(in) :: r(3),q(3,3),rion
    real(rk), intent(out) :: foot(3),df(3,3)
    integer,intent(out) :: hemisphere,status
    real(rk) :: rm(3),fm(3),d(3,3),rr,r2,rho2,a,g,dg(3),sgn
    integer :: i,j
    foot=0;df=0;hemisphere=0;status=1
    if(.not.all(ieee_is_finite(r)).or..not.all(ieee_is_finite(q))) return
    if(.not.ieee_is_finite(rion).or.rion<=0) return
    rm=matmul(transpose(q),r);r2=sum(rm*rm);rr=sqrt(r2)
    if(rr<=rion) return
    if(abs(rm(3))<=32*epsilon(rr)*rr) then
      status=2;return
    endif
    hemisphere=merge(1,-1,rm(3)>0)
    sgn=real(hemisphere,rk);rho2=sum(rm(:2)**2)
    a=(rion/rr)**1.5_rk
    g=1-rion*rho2/rr**3
    if(g<=0) return
    fm=[a*rm(1),a*rm(2),sgn*rion*sqrt(g)]
    d=0
    do i=1,2
      do j=1,3
        d(i,j)=-1.5_rk*a*rm(i)*rm(j)/r2
      enddo
      d(i,i)=d(i,i)+a
    enddo
    dg=3*rion*rho2*rm/rr**5
    dg(:2)=dg(:2)-2*rion*rm(:2)/rr**3
    d(3,:)=sgn*rion*dg/(2*sqrt(g))
    foot=matmul(q,fm)
    df=matmul(q,matmul(d,transpose(q)))
    status=0
  end subroutine

  pure subroutine mi_outward_point(foot,q,rion,rsample,r,status)
    ! Analytic inverse to a SPHERICAL sampling shell. Caller must verify the
    ! actual MHD shell geometry; never substitute its mean radius silently.
    real(rk), intent(in) :: foot(3),q(3,3),rion,rsample
    real(rk), intent(out) :: r(3)
    integer,intent(out) :: status
    real(rk)::fm(3),rm(3),rho2
    r=0;status=1
    if(rion<=0.or.rsample<=rion) return
    if(abs(norm2(foot)/rion-1)>1.e-10_rk) return
    fm=matmul(transpose(q),foot)
    rm(:2)=fm(:2)*(rsample/rion)**1.5_rk
    rho2=sum(rm(:2)**2)
    if(rho2>rsample**2*(1+1.e-12_rk)) return
    rm(3)=sign(sqrt(max(0._rk,rsample**2-rho2)),fm(3))
    r=matmul(q,rm);status=0
  end subroutine

  pure subroutine mi_fac_radial(j_source,b_source,b_foot,foot,jr,status)
    ! J_source [A/m^2], B at source AND footprint [T]. No constant polar-B
    ! approximation, no magnitude-only sign convention. Assumes j_parallel/B
    ! conserved over this inner dipolar segment; does not trace the outer MHD.
    real(rk),intent(in)::j_source(3),b_source(3),b_foot(3),foot(3)
    real(rk),intent(out)::jr
    integer,intent(out)::status
    real(rk)::bsq,rf
    jr=0;status=1;bsq=sum(b_source*b_source);rf=norm2(foot)
    if(bsq<=tiny(bsq).or.rf<=0) return
    jr=dot_product(j_source,b_source)/bsq*dot_product(b_foot,foot/rf)
    if(.not.ieee_is_finite(jr)) return
    status=0
  end subroutine

  pure function mi_pullback_e(df,gradient_foot) result(e)
    ! grad Phi at foot in V/m, resulting physical E at source in V/m.
    ! Rotating a previously computed E alone would omit the mapping Jacobian.
    real(rk),intent(in)::df(3,3),gradient_foot(3)
    real(rk)::e(3)
    e=-matmul(transpose(df),gradient_foot)
  end function

  pure function mi_edge_lfm(phi_a,phi_b) result(emf)
    ! LFM evolves a +curl (-c E) edge integral. Potential inputs in volts;
    ! this is the donor 1e8*(Phi_b-Phi_a), NOT physical integral E.dl.
    real(rk),intent(in)::phi_a,phi_b
    real(rk)::emf
    emf=1.e8_rk*(phi_b-phi_a)
  end function

  pure subroutine mi_wall_drift(e_t,b,n,v_normal,bn_relative_tolerance,v,status)
    ! Ideal boundary velocity constrained by a specified normal velocity:
    ! E_t = B_n n x v_t - v_n n x B_t.
    ! No arbitrary projection of ExB: such projection changes E_t when B
    ! is tilted. Near B_n=0, return status=2 and require a verified analytic
    ! cap-rim limit; never hide this singularity with a velocity/B floor.
    real(rk),intent(in)::e_t(3),b(3),n(3),v_normal,bn_relative_tolerance
    real(rk),intent(out)::v(3)
    integer,intent(out)::status
    real(rk)::bn,bmag,bt(3),et(3)
    v=0;status=1;bmag=norm2(b)
    if(abs(norm2(n)-1)>1.e-10_rk.or.bmag<=0.or.bn_relative_tolerance<0) return
    bn=dot_product(b,n);bt=b-bn*n
    et=e_t-dot_product(e_t,n)*n
    if(abs(bn)<=bn_relative_tolerance*bmag) then
      status=2;return
    endif
    v=(v_normal*bt-cross3(n,et))/bn+v_normal*n
    if(.not.all(ieee_is_finite(v))) return
    status=0
  end subroutine
end module
