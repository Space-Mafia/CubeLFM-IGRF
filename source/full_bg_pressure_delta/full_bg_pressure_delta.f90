subroutine diagnostic_full_pressure_delta(v,face_area_vector,delta_flux)
 use iso_fortran_env,only:real64
 use lfm_cube_mesh,only:cube_cross
 implicit none
 ! Exact previously qualified geometry correction; only static B0 callbacks bind to IGRF.
  ! get_face orders its corners around the perimeter: 00,10,11,01.
  ! Return the area-integrated magnetic-momentum flux correction in cgs.
  ! Original G2INT quadrature nodes/weights and B0 functions are retained.
  real(real64),intent(in)::v(3,4),face_area_vector(3)
  real(real64),intent(out)::delta_flux(3)
  real(real64)::normal(3),de(3),dp(3),dep(3),point(3),da(3),b0(3)
  real(real64)::eta,psi,weight,pressure,pi4inv
  real(real64),external::diagnostic_full_bx,diagnostic_full_by,diagnostic_full_bz
  real(real64),parameter::node(12)=[ &
   .1252334085_real64,.3678314989_real64,.5873179542_real64, &
   .7699026741_real64,.9041172563_real64,.9815606342_real64, &
   -.1252334085_real64,-.3678314989_real64,-.5873179542_real64, &
   -.7699026741_real64,-.9041172563_real64,-.9815606342_real64]
  real(real64),parameter::wt(12)=[ &
   .2491470458_real64,.2334925365_real64,.2031674267_real64, &
   .1600783285_real64,.1069393259_real64,.0471753363_real64, &
   .2491470458_real64,.2334925365_real64,.2031674267_real64, &
   .1600783285_real64,.1069393259_real64,.0471753363_real64]
  integer::ig,jg
  normal=face_area_vector/norm2(face_area_vector)
  de=v(:,2)-v(:,1);dp=v(:,4)-v(:,1);dep=v(:,3)+v(:,1)-v(:,2)-v(:,4)
  delta_flux=0;pi4inv=1/(4*acos(-1._real64))
  do ig=1,12
   eta=(1+node(ig))/2
   do jg=1,12
    psi=(1+node(jg))/2;weight=wt(ig)*wt(jg)/4
    point=v(:,1)+eta*de+psi*dp+eta*psi*dep
    da=cube_cross(de+psi*dep,dp+eta*dep)
    b0=[diagnostic_full_bx(point(1),point(2),point(3)), &
     diagnostic_full_by(point(1),point(2),point(3)),diagnostic_full_bz(point(1),point(2),point(3))]
    pressure=dot_product(b0,b0)/2
    ! Correct only the B0-pressure normal geometry. Original magnetic tension,
    ! residual/cross-field stress, kinetic diffusion and Boris algebra remain.
    delta_flux=delta_flux+weight*pi4inv*pressure*(da-normal*norm2(da))
   enddo
  enddo
end subroutine
