! Units/context adapter only; tagged C field/A and native GINT/G2INT algorithms.
module lfm_cube_yy_background
 use iso_c_binding,only:c_double,c_int
 use lfm_cube_yy_igrf,only:cube_yy_igrf_configure
 implicit none
 private
 public::cube_yy_background_configure,cube_yy_background_context,cube_yy_background_ready
 public::cube_yy_evaluate_cgs,cube_yy_context_potential,cube_yy_static_background
 public::cube_yy_ct_evaluate
 real(c_double),save::length_cm=0,field_gauss=0,time_unit=0,context_time=0
 integer(c_int),save::background_role=0,context_role=0
 logical,save::ready=.false.
 interface
  integer(c_int) function cube_yy_ct_configure(length_m,magnetic_tesla,time_s,omega,phase) &
    bind(C,name='cube_yy_ct_configure')
   import c_double,c_int
   real(c_double),value::length_m,magnetic_tesla,time_s,omega,phase
  end function
  integer(c_int) function cube_yy_ct_evaluate(xyz,time,role,taper,b,a) bind(C,name='cube_yy_ct_evaluate')
   import c_double,c_int
   real(c_double),intent(in)::xyz(3)
   real(c_double),value::time
   integer(c_int),value::role,taper
   real(c_double),intent(out)::b(3),a(3)
  end function
 end interface
contains
 subroutine cube_yy_background_configure(length_m,magnetic_tesla,time_s,omega,phase,status)
  real(c_double),intent(in)::length_m,magnetic_tesla,time_s,omega,phase
  integer,intent(out)::status
  status=cube_yy_ct_configure(length_m,magnetic_tesla,time_s,omega,phase)
  if(status/=0) return
  ! Reuse the already qualified uncut field/tracer object with identical inputs.
  status=cube_yy_igrf_configure(length_m,magnetic_tesla,omega,phase)
  if(status/=0) return
  length_cm=100*length_m;field_gauss=1.e4_c_double*magnetic_tesla;time_unit=time_s
  ! This is exactly the tagged background callback's role selection.
  background_role=0;if(omega>0) background_role=1
  context_time=0;context_role=0;ready=.true.
 end subroutine
 logical function cube_yy_background_ready()
  cube_yy_background_ready=ready
 end function
 subroutine cube_yy_background_context(time_s,role)
  real(c_double),intent(in)::time_s
  integer,intent(in)::role
  if(.not.ready.or.(role/=0.and.role/=1)) error stop 'YY background context not configured'
  context_time=time_s
  context_role=0;if(role==1) context_role=background_role
 end subroutine
 subroutine cube_yy_evaluate_cgs(x,time_s,role,taper,b,a,status)
  real(c_double),intent(in)::x(3),time_s
  integer,intent(in)::role,taper
  real(c_double),intent(out)::b(3),a(3)
  integer,intent(out)::status
  real(c_double)::bc(3),ac(3)
  if(.not.ready) error stop 'YY background units not configured'
  status=cube_yy_ct_evaluate(x/length_cm,time_s/time_unit,int(role,c_int),int(taper,c_int),bc,ac)
  if(status/=0) return
  b=bc*field_gauss;a=ac*field_gauss*length_cm
 end subroutine
 subroutine cube_yy_context_potential(x,a)
  real(c_double),intent(in)::x(3)
  real(c_double),intent(out)::a(3)
  real(c_double)::b(3)
  integer::st
  call cube_yy_evaluate_cgs(x,context_time,int(context_role),1,b,a,st)
  if(st/=0) error stop 'Tagged YY A callback failed'
 end subroutine
 subroutine cube_yy_static_background(x,b)
  real(c_double),intent(in)::x(3)
  real(c_double),intent(out)::b(3)
  real(c_double)::a(3)
  integer::st
  call cube_yy_evaluate_cgs(x,0._c_double,int(background_role),1,b,a,st)
  if(st/=0) error stop 'Tagged YY B0 callback failed'
 end subroutine
end module

! Legacy scalar callback ABI retained for existing native GINT/G2INT.
function cube_yy_ax(x,y,z) result(v)
 use iso_c_binding,only:c_double
 use lfm_cube_yy_background,only:cube_yy_context_potential
 real(c_double)::x,y,z,v,a(3)
 call cube_yy_context_potential([x,y,z],a);v=a(1)
end function
function cube_yy_ay(x,y,z) result(v)
 use iso_c_binding,only:c_double
 use lfm_cube_yy_background,only:cube_yy_context_potential
 real(c_double)::x,y,z,v,a(3)
 call cube_yy_context_potential([x,y,z],a);v=a(2)
end function
function cube_yy_az(x,y,z) result(v)
 use iso_c_binding,only:c_double
 use lfm_cube_yy_background,only:cube_yy_context_potential
 real(c_double)::x,y,z,v,a(3)
 call cube_yy_context_potential([x,y,z],a);v=a(3)
end function
function cube_yy_bxqq0(x,y,z) result(v)
 use iso_c_binding,only:c_double
 use lfm_cube_yy_background,only:cube_yy_static_background
 real(c_double)::x,y,z,v,b(3)
 call cube_yy_static_background([x,y,z],b);v=b(1)
end function
function cube_yy_byqq0(x,y,z) result(v)
 use iso_c_binding,only:c_double
 use lfm_cube_yy_background,only:cube_yy_static_background
 real(c_double)::x,y,z,v,b(3)
 call cube_yy_static_background([x,y,z],b);v=b(2)
end function
function cube_yy_bzqq0(x,y,z) result(v)
 use iso_c_binding,only:c_double
 use lfm_cube_yy_background,only:cube_yy_static_background
 real(c_double)::x,y,z,v,b(3)
 call cube_yy_static_background([x,y,z],b);v=b(3)
end function
