module lfm_cube_yy_native_cap_frame
 use iso_fortran_env,only:real64
 use lfm_rotating_mi_geometry,only:mi_frame
 implicit none
 private
 public::yy_native_cap_frame
 include 'yy_rotation_profile.inc'
contains
 function yy_native_cap_frame(time)result(q)
  real(real64),intent(in)::time
  real(real64)::q(3,3),identity(3,3),axis(3),alpha,phase
  integer::i
  ! Magnetic north axis from the SAME frozen IGRF degree1 coefficients.
  ! Reuse LFM mi_frame; no change to the full degree5 YY field/tracer.
  axis=-[-1450.9_real64,4652.5_real64,-29404.8_real64];axis=axis/norm2(axis)
  alpha=acos(axis(3));phase=atan2(axis(2),axis(1))
  identity=0;do i=1,3;identity(i,i)=1;enddo
  q=mi_frame(identity,alpha,phase+yy_spin_phase0_rad,yy_spin_omega_rad_s,time)
 end function
end module
