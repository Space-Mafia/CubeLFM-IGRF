! Small magnetic-frame boundary state, identical on every I=0 MPI tile.
! No MHD field/corner arrays are stored here; no Cartesian E is held in time.
module lfm_rotating_mi_state
  use iso_fortran_env,only:real64
  use, intrinsic::ieee_arithmetic,only:ieee_is_finite
  use lfm_tilt_geometry,only:tilt_frame,tilt_alpha,tilt_phase0,tilt_omega,tilt_mu,dipole_field
  use lfm_rotating_mi_geometry,only:mi_frame,mi_wall_drift
  use lfm_rotating_mi_map,only:mi_inner_potential
  implicit none
  private
  integer,parameter::rk=real64
  logical,public,save::rotating_mi=.false.,rotating_mi_zero_source=.false.
  real(rk),allocatable,save::angles(:),azimuths(:),potential(:,:)
  real(rk),save::radius=0,held_time=0
  public::mi_state_install,mi_state_export,mi_potential_point,mi_velocity_point
  public::mi_state_validate,mi_state_time
contains
  subroutine mi_state_validate(buffer,n,m,status)
    integer,intent(in)::n,m
    real(rk),intent(in)::buffer(:)
    integer,intent(out)::status
    real(rk)::tab(n,m),th(n),ph(m)
    status=1
    if(n<3.or.m<5.or.mod(m-1,2)/=0.or.size(buffer)/=2+n+m+n*m) return
    if(.not.all(ieee_is_finite(buffer)).or.buffer(1)<=0.or.buffer(2)<0) return
    th=buffer(3:2+n);ph=buffer(3+n:2+n+m)
    tab=reshape(buffer(3+n+m:),[n,m])
    if(any(th(2:)<=th(:n-1)).or.any(ph(2:)<=ph(:m-1))) return
    if(abs(ph(m)-ph(1)-8*atan(1._rk))>1.e-10_rk) return
    if(any(tab(:,1)/=0).or.any(tab(:,(m+1)/2)/=0).or.any(tab(:,m)/=0)) return
    if(any(tab(1,:)/=0).or.any(tab(n,:)/=0)) return
    status=0
  end subroutine
  subroutine mi_state_install(buffer,n,m,status)
    integer,intent(in)::n,m
    real(rk),intent(in)::buffer(:)
    integer,intent(out)::status
    call mi_state_validate(buffer,n,m,status)
    if(status/=0) return
    angles=buffer(3:2+n);azimuths=buffer(3+n:2+n+m)
    potential=reshape(buffer(3+n+m:),[n,m])
    radius=buffer(1);held_time=buffer(2)
  end subroutine
  function mi_state_time() result(t)
    real(rk)::t
    t=held_time
  end function
  subroutine mi_state_export(buffer,status)
    real(rk),intent(out)::buffer(:)
    integer,intent(out)::status
    integer::n,m
    status=1
    if(.not.allocated(potential)) return
    n=size(angles);m=size(azimuths)
    if(size(buffer)/=2+n+m+n*m) return
    buffer=[radius,held_time,angles,azimuths,reshape(potential,[n*m])]
    status=0
  end subroutine
  subroutine point_e(r_cm,time,value,e,q,status)
    real(rk),intent(in)::r_cm(3),time
    real(rk),intent(out)::value,e(3),q(3,3)
    integer,intent(out)::status
    value=0;e=0;status=0
    q=mi_frame(tilt_frame,tilt_alpha,tilt_phase0,tilt_omega,time)
    if(.not.rotating_mi.or..not.allocated(potential)) return
    call mi_inner_potential(.01_rk*r_cm,q,radius,angles,azimuths,potential,value,e,status)
    ! Below-wall closed lines have no magnetospheric extension. The installed
    ! complete table has exact equator/axis zeros and regular endpoint slopes.
    if(status==2.and.norm2(.01_rk*r_cm)<=radius*(1+1.e-10_rk)) then
      value=0;e=0;status=0
    endif
  end subroutine
  subroutine mi_potential_point(r_cm,time,value,status)
    real(rk),intent(in)::r_cm(3),time
    real(rk),intent(out)::value
    integer,intent(out)::status
    real(rk)::e(3),q(3,3)
    call point_e(r_cm,time,value,e,q,status)
  end subroutine
  subroutine mi_velocity_point(r_cm,time,v_cm,status)
    real(rk),intent(in)::r_cm(3),time
    real(rk),intent(out)::v_cm(3)
    integer,intent(out)::status
    real(rk)::value,e(3),q(3,3),b(3),n(3)
    v_cm=0
    call point_e(r_cm,time,value,e,q,status)
    if(status/=0.or.all(e==0)) return
    b=dipole_field(.01_rk*r_cm,tilt_mu*1.e-10_rk*q(:,3));n=r_cm/norm2(r_cm)
    call mi_wall_drift(e,b,n,0._rk,1.e-14_rk,v_cm,status)
    v_cm=100*v_cm
  end subroutine
end module
