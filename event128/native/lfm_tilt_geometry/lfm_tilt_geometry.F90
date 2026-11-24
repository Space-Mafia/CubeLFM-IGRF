! Coordinate-independent tilted dipole tools, plus opt-in driver configuration.
! Cartesian r and moment use matching length units; moment is B*length**3.
! LFM convention: moment magnitude M = -GEOQMU for the current +Z Jupiter.
module lfm_tilt_geometry
  use iso_fortran_env, only: real64
  use iso_c_binding, only: c_double
  implicit none
  private
  public :: dipole_rotation, dipole_moment, dipole_field, dipole_potential
  public :: dipole_footpoint, edge_potential, rotation_edge_emf, cross3
  public :: tilt_edge_c
  public :: tilt_configure, tilt_enabled, tilt_mu, tilt_alpha, tilt_omega, tilt_phase0
  public :: tilt_frame
  real(real64), save :: tilt_frame(3,3)=reshape([1._real64,0._real64,0._real64, &
    0._real64,1._real64,0._real64,0._real64,0._real64,1._real64],[3,3])
  logical, save :: tilt_enabled=.false.
  real(real64), save :: tilt_mu=0,tilt_alpha=0,tilt_omega=0,tilt_phase0=0
  real(real64), parameter :: pi=acos(-1._real64)
contains
  subroutine tilt_configure(radius_cm,beq_g,angle_deg,period_s,phase_deg,spin_deg,spin_az_deg,boundary_rotation)
    real(real64),intent(in)::radius_cm,beq_g,angle_deg,period_s,phase_deg
    real(real64),intent(in),optional::spin_deg,spin_az_deg
    logical,intent(in),optional::boundary_rotation
    real(real64)::beta,az
    beta=0;az=0
    if(present(spin_deg)) beta=spin_deg*pi/180
    if(present(spin_az_deg)) az=spin_az_deg*pi/180
    tilt_frame=dipole_rotation(beta,az)
    tilt_enabled=abs(angle_deg)>0
    tilt_enabled=tilt_enabled.or.abs(beta)>0
    ! Explicit rigid-boundary mode must also work for both inclinations zero.
    if(present(boundary_rotation)) tilt_enabled=tilt_enabled.or.boundary_rotation
    tilt_mu=-beq_g*radius_cm**3
    tilt_alpha=angle_deg*pi/180
    tilt_phase0=phase_deg*pi/180
    tilt_omega=0
    if(period_s>0) tilt_omega=2*pi/period_s
  end subroutine
  pure function cross3(a,b) result(c)
    real(real64), intent(in) :: a(3),b(3)
    real(real64) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  pure function dipole_rotation(alpha,phase) result(q)
    ! Magnetic coordinates -> inertial coordinates: Rz(phase)*Ry(alpha).
    ! alpha>0 and phase=0 tilts the +Z magnetic axis toward +X.
    real(real64), intent(in) :: alpha,phase
    real(real64) :: q(3,3),ca,sa,cp,sp
    ca=cos(alpha);sa=sin(alpha);cp=cos(phase);sp=sin(phase)
    q(:,1)=[cp*ca,sp*ca,-sa]
    q(:,2)=[-sp,cp,0._real64]
    q(:,3)=[cp*sa,sp*sa,ca]
  end function

  pure function dipole_moment(mu,alpha,omega,phase0,time,frame) result(m)
    real(real64), intent(in) :: mu,alpha,omega,phase0,time
    real(real64), intent(in),optional :: frame(3,3)
    real(real64) :: m(3),q(3,3)
    q=dipole_rotation(alpha,modulo(phase0+omega*time,2*pi))
    m=mu*q(:,3)
    if(present(frame)) m=matmul(frame,m)
  end function

  pure function dipole_field(r,m) result(b)
    real(real64), intent(in) :: r(3),m(3)
    real(real64) :: b(3),r2
    r2=sum(r*r)
    b=(3*r*dot_product(m,r)-m*r2)/(r2*r2*sqrt(r2))
  end function

  pure function dipole_potential(r,m) result(a)
    real(real64), intent(in) :: r(3),m(3)
    real(real64) :: a(3),r2
    r2=sum(r*r)
    a=cross3(m,r)/(r2*sqrt(r2))
  end function

  pure subroutine dipole_footpoint(r,mu,alpha,phase,rion,foot,b_ratio,status)
    ! Pure dipole inward mapping. The magnetic hemisphere follows r dot mhat.
    ! b_ratio=|B_source|/|B_foot| (FAC mapping must choose its own direction).
    real(real64), intent(in) :: r(3),mu,alpha,phase,rion
    real(real64), intent(out) :: foot(3),b_ratio
    integer, intent(out) :: status
    real(real64) :: q(3,3),rm(3),fm(3),rr,rho,sion,cphi,sphi,m(3)
    q=dipole_rotation(alpha,phase)
    rm=matmul(transpose(q),r);rr=norm2(rm)
    foot=0;b_ratio=0;status=1
    if(rion<=0.or.rr<rion.or.mu==0) return
    rho=norm2(rm(:2))
    if(rho>epsilon(rr)*rr) then
      cphi=rm(1)/rho;sphi=rm(2)/rho
    else
      cphi=1;sphi=0
    endif
    sion=min(1._real64,sqrt(rion/rr)*rho/rr)
    fm=rion*[sion*cphi,sion*sphi,sign(sqrt(max(0._real64,1-sion*sion)),rm(3))]
    foot=matmul(q,fm)
    m=mu*q(:,3)
    b_ratio=norm2(dipole_field(r,m))/norm2(dipole_field(foot,m))
    status=0
  end subroutine

  pure function edge_potential(a,b,m) result(value)
    ! Straight physical corner edge, 8-point Gauss-Legendre integral A.dl.
    ! Edges must not pass through r=0; physical LFM inner boundary excludes it.
    real(real64), intent(in) :: a(3),b(3),m(3)
    real(real64) :: value,r(3),d(3)
    real(real64), parameter :: nodes(8)=[ &
      -.9602898564975363_real64,-.7966664774136267_real64, &
      -.5255324099163290_real64,-.1834346424956498_real64, &
       .1834346424956498_real64, .5255324099163290_real64, &
       .7966664774136267_real64, .9602898564975363_real64]
    real(real64), parameter :: weights(8)=[ &
      .1012285362903763_real64,.2223810344533745_real64, &
      .3137066458778873_real64,.3626837833783620_real64, &
      .3626837833783620_real64,.3137066458778873_real64, &
      .2223810344533745_real64,.1012285362903763_real64]
    integer :: n
    d=b-a;value=0
    do n=1,8
      r=.5_real64*(a+b+nodes(n)*d)
      value=value+.5_real64*weights(n)*dot_product(dipole_potential(r,m),d)
    enddo
  end function

  pure function sinc(x) result(y)
    real(real64), intent(in) :: x
    real(real64) :: y
    if(abs(x)<1.e-4_real64) then
      y=1-x*x/6+x**4/120-x**6/5040
    else
      y=sin(x)/x
    endif
  end function

  pure function rotation_edge_emf(a,b,mu,alpha,omega,phase0,time,dt,frame) result(u)
    ! Positive LFM curl variable, NOT physical SI E.
    ! For rigid rotation: u x B = dA/dt + grad(u dot A).
    ! Integrate the A difference over the edge and the scalar gradient exactly.
    ! Time averages are analytic; CT curl of endpoint scalar differences cancels.
    ! No dt alignment, no subtraction of nearly equal large edge potentials.
    real(real64), intent(in) :: a(3),b(3),mu,alpha,omega,phase0,time,dt
    real(real64), intent(in),optional :: frame(3,3)
    real(real64) :: u,mavg(3),mdotavg(3),mid(3),w(3),h,psia,psib
    h=.5_real64*omega*dt
    mid=dipole_moment(mu,alpha,omega,phase0,time+.5_real64*dt)
    mavg=mid;mavg(:2)=mid(:2)*sinc(h)
    w=[0._real64,0._real64,omega]
    ! Average in the spin frame, then rotate ALL vectors into Cartesian axes.
    if(present(frame)) then
      mavg=matmul(frame,mavg)
      w=matmul(frame,w)
    endif
    mdotavg=cross3(w,mavg)
    psia=dot_product(cross3(w,a),dipole_potential(a,mavg))
    psib=dot_product(cross3(w,b),dipole_potential(b,mavg))
    u=edge_potential(a,b,mdotavg)+psib-psia
  end function

  subroutine tilt_edge_c(a,b,mu,alpha,omega,phase0,time,dt,ai0,ai1,u) bind(C)
    ! Small testing ABI; callers supply double, lengths in RJ are fine for tests.
    real(c_double), intent(in) :: a(3),b(3)
    real(c_double), value :: mu,alpha,omega,phase0,time,dt
    real(c_double), intent(out) :: ai0,ai1,u
    ai0=edge_potential(a,b,dipole_moment(mu,alpha,omega,phase0,time))
    ai1=edge_potential(a,b,dipole_moment(mu,alpha,omega,phase0,time+dt))
    u=rotation_edge_emf(a,b,mu,alpha,omega,phase0,time,dt)
  end subroutine
end module
