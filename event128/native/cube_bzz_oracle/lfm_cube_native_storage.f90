! Rank-local storage ABI for the unchanged native metric/BZZ/CFL/CT routines.
! No field-recovery arithmetic is replaced in this oracle.
module lfm_cube_native_storage
  use iso_fortran_env,only:real64
  implicit none
  integer :: nl_i,nl_j,nl_k,nl_ip1,nl_jp1,nl_kp1
  integer,parameter :: ilow=1,ihigh=65,jlow=1,jhigh=33,klow=1,khigh=33
  logical :: mhd_pe=.true.
  real(real64),save :: x(65,33,33),y(65,33,33),z(65,33,33)
  real(real64),save :: face_i(65,33,33),face_j(65,33,33),face_k(65,33,33),volume(65,33,33)
  real(real64),save :: volinv(65,33,33),cell_length(65,33,33)
  real(real64),save :: edge_i(65,33,33),edge_j(65,33,33),edge_k(65,33,33)
  real(real64),save :: tran_i(65,33,33,9),tran_j(65,33,33,9),tran_k(65,33,33,9)
  real(real64),save :: vec_i(65,33,33,9),vec_j(65,33,33,9),vec_k(65,33,33,9)
  real(real64),save :: rho_g(65,33,33),c_g(65,33,33),vx_g(65,33,33),vy_g(65,33,33),vz_g(65,33,33)
  real(real64),save :: bx_g(65,33,33),by_g(65,33,33),bz_g(65,33,33),bsqqk(65,33,33)
  real(real64) :: pi4inv,cainv
  integer :: i_global_off=0,j_global_off=0,k_global_off=0
  real(real64),save :: bi_g(65,33,33),bj_g(65,33,33),bk_g(65,33,33)
  real(real64),save :: bi_dbl(65,33,33),bj_dbl(65,33,33),bk_dbl(65,33,33)
  real(real64),save :: ei_g(65,33,33),ej_g(65,33,33),ek_g(65,33,33)
end module
