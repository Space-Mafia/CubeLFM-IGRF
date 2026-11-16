module diagnostic_full_A_bg
 use iso_c_binding,only:c_double
 use iso_fortran_env,only:real64
 use lfm_cube_yy_background,only:cube_yy_evaluate_cgs
 implicit none
 real(real64),save::diagnostic_epoch=0
 interface
 subroutine memo(x,t,b)bind(C,name='cube_yy_hydro_potential_memo')
 import c_double
 real(c_double),intent(in)::x(3)
 real(c_double),value::t
 real(c_double),intent(out)::b(3)
 end subroutine
 subroutine reset()bind(C,name='cube_yy_hydro_potential_reset')
 end subroutine
 end interface
 contains
 function component(x,y,z,i)result(v)
 real(real64),intent(in)::x,y,z
 integer,intent(in)::i
 real(real64)::v,b(3)
 call memo([x,y,z],diagnostic_epoch,b);v=b(i)
 end function
end module
function diagnostic_full_ax(x,y,z)result(v)
 use diagnostic_full_A_bg
 real(real64)::x,y,z,v
 v=component(x,y,z,1)
end function
function diagnostic_full_ay(x,y,z)result(v)
 use diagnostic_full_A_bg
 real(real64)::x,y,z,v
 v=component(x,y,z,2)
end function
function diagnostic_full_az(x,y,z)result(v)
 use diagnostic_full_A_bg
 real(real64)::x,y,z,v
 v=component(x,y,z,3)
end function
subroutine diagnostic_full_A_set_time(t)
 use diagnostic_full_A_bg
 real(real64),intent(in)::t
 diagnostic_epoch=t;call reset()
end subroutine
subroutine full_potential_uncached(x,t,b)bind(C,name='cube_yy_hydro_potential_uncached')
 use iso_c_binding,only:c_double
 use lfm_cube_yy_background,only:cube_yy_evaluate_cgs
 real(c_double),intent(in)::x(3),t
 real(c_double),intent(out)::b(3)
 real(c_double)::a(3)
 integer::st
 call cube_yy_evaluate_cgs(x,t,0,1,a,b,st)
 if(st/=0)error stop 'existing full prescribed field callback'
end subroutine
