! C ABI only: no new field synthesis, line tracer, solver or rotation algorithm.
module lfm_cube_yy_igrf
 use iso_c_binding, only:c_double,c_int
 implicit none
 private
 public::cube_yy_igrf_configure,cube_yy_igrf_field,cube_yy_igrf_trace
 interface
  integer(c_int) function cube_yy_igrf_configure(length_m,magnetic_tesla,omega,phase) &
    bind(C,name='cube_yy_igrf_configure')
   import c_double,c_int
   real(c_double),value::length_m,magnetic_tesla,omega,phase
  end function
  integer(c_int) function cube_yy_igrf_field(xyz,time_s,direction,magnitude) &
    bind(C,name='cube_yy_igrf_field')
   import c_double,c_int
   real(c_double),intent(in)::xyz(3)
   real(c_double),value::time_s
   real(c_double),intent(out)::direction(3),magnitude
  end function
  integer(c_int) function cube_yy_igrf_trace(spherical_start,target_radius,time_s,endpoint) &
    bind(C,name='cube_yy_igrf_trace')
   import c_double,c_int
   real(c_double),intent(in)::spherical_start(3)
   real(c_double),value::target_radius,time_s
   real(c_double),intent(out)::endpoint(4)
  end function
 end interface
end module
