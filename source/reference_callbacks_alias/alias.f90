! ABI forwarding only; underlying ct_phys provider qualified on actual Cube faces.
subroutine diagnostic_full_set_time(t)
 use iso_fortran_env,only:real64
 real(real64),intent(in)::t
 call ct_phys_set_time(t)
end subroutine
function diagnostic_full_bx(x,y,z) result(v)
 use iso_fortran_env,only:real64
 real(real64)::x,y,z,v
 real(real64),external::ct_phys_bx
 v=ct_phys_bx(x,y,z)
end function
function diagnostic_full_by(x,y,z) result(v)
 use iso_fortran_env,only:real64
 real(real64)::x,y,z,v
 real(real64),external::ct_phys_by
 v=ct_phys_by(x,y,z)
end function
function diagnostic_full_bz(x,y,z) result(v)
 use iso_fortran_env,only:real64
 real(real64)::x,y,z,v
 real(real64),external::ct_phys_bz
 v=ct_phys_bz(x,y,z)
end function
