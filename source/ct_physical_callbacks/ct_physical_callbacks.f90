module ct_physical_bg
 use iso_c_binding,only:c_double
 use iso_fortran_env,only:real64
 use lfm_cube_yy_background,only:cube_yy_evaluate_cgs
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 real(real64),save::epoch=0,length_cm=0
 interface
 function getter()bind(C,name='cube_yy_ctphys_length_cm')result(v)
 import c_double
 real(c_double)::v
 end function
 subroutine weight(r,w,dw)bind(C,name='cube_yy_ctphys_taper')
 import c_double
 real(c_double),value::r
 real(c_double),intent(out)::w,dw
 end subroutine
 subroutine memo(x,t,b)bind(C,name='cube_yy_ctphys_vector_memo')
 import c_double
 real(c_double),intent(in)::x(3)
 real(c_double),value::t
 real(c_double),intent(out)::b(3)
 end subroutine
 subroutine reset()bind(C,name='cube_yy_ctphys_vector_reset')
 end subroutine
 end interface
 contains
 function component(x,y,z,i)result(v)
 real(real64),intent(in)::x,y,z
 integer,intent(in)::i
 real(real64)::v,b(3)
 call memo([x,y,z],epoch,b);v=b(i)
 end function
end module
subroutine ct_phys_set_time(t)
 use ct_physical_bg
 real(real64),intent(in)::t
 length_cm=getter();if(.not.ieee_is_finite(length_cm).or.length_cm<=0)error stop 'original CT normalization length'
 epoch=t;call reset()
end subroutine
function ct_phys_bx(x,y,z)result(v)
 use ct_physical_bg
 real(real64)::x,y,z,v
 v=component(x,y,z,1)
end function
function ct_phys_by(x,y,z)result(v)
 use ct_physical_bg
 real(real64)::x,y,z,v
 v=component(x,y,z,2)
end function
function ct_phys_bz(x,y,z)result(v)
 use ct_physical_bg
 real(real64)::x,y,z,v
 v=component(x,y,z,3)
end function
subroutine ct_physical_uncached(x,t,b)bind(C,name='cube_yy_ctphys_vector_uncached')
 use ct_physical_bg,only:c_double,real64,length_cm,weight,cube_yy_evaluate_cgs
 real(c_double),intent(in)::x(3),t
 real(c_double),intent(out)::b(3)
 real(c_double)::a(3),r,w,dw,gradient(3),cross(3)
 integer::st
 r=norm2(x)
 if(r<=0.or.length_cm<=0)error stop 'original CT field point/context'
 call weight(r/length_cm,w,dw)
 if(dw==0)then
  ! Exactly retain original released provider in constant-weight regions.
  call cube_yy_evaluate_cgs(x,t,0,1,b,a,st)
 else
  call cube_yy_evaluate_cgs(x,t,0,0,b,a,st)
  gradient=(dw/length_cm)*(x/r)
  cross=[gradient(2)*a(3)-gradient(3)*a(2),gradient(3)*a(1)-gradient(1)*a(3),gradient(1)*a(2)-gradient(2)*a(1)]
  b=w*b+cross
 endif
 if(st/=0)error stop 'original full field/A provider'
end subroutine
