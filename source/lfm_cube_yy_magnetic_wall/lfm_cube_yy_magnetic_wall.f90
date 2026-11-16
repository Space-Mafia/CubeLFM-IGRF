! Derived geometry only; no new checkpoint state or model-type member.
module lfm_cube_yy_magnetic_wall
 use iso_c_binding,only:c_ptr,c_null_ptr,c_associated,c_int,c_double
 implicit none
 private
 public::yy_magnetic_wall_prepare,yy_magnetic_wall_apply
 type(c_ptr),save::wall=c_null_ptr
 interface
  integer(c_int) function c_prepare(nj,nk,g,vertices,area,norms,handle) &
    bind(C,name='cube_yy_magnetic_wall_create')
   import c_int,c_double,c_ptr
   integer(c_int),value::nj,nk,g
   real(c_double),intent(in)::vertices(*),area(*),norms(*)
   type(c_ptr),intent(out)::handle
  end function
  integer(c_int) function c_apply(handle,ni,time,dt_old,tfac,bf) &
    bind(C,name='cube_yy_magnetic_wall_apply')
   import c_int,c_double,c_ptr
   type(c_ptr),value::handle
   integer(c_int),value::ni
   real(c_double),value::time,dt_old,tfac
   real(c_double),intent(inout)::bf(*)
  end function
  subroutine c_destroy(handle)bind(C,name='cube_yy_magnetic_wall_destroy')
   import c_ptr
   type(c_ptr),value::handle
  end subroutine
 end interface
contains
 subroutine yy_magnetic_wall_prepare(g,n,vertices,area,norms,status)
  integer,intent(in)::g,n(3)
  real(c_double),intent(in)::vertices(:,-g:,-g:,-g:)
  real(c_double),intent(in)::area(:,1-g:,1-g:,1-g:,:),norms(1-g:,1-g:,1-g:,:)
  integer,intent(out)::status
  real(c_double),allocatable::v(:,:,:,:),a(:,:,:,:,:),s(:,:,:,:)
  if(c_associated(wall))call c_destroy(wall)
  wall=c_null_ptr
  ! One initialization copy. The hot apply passes the full contiguous bf.
  v=vertices(:,-g:1,-g:n(2)+g,-g:n(3)+g)
  a=area(:,1-g:1,1-g:n(2)+g,1-g:n(3)+g,:)
  s=norms(1-g:1,1-g:n(2)+g,1-g:n(3)+g,:)
  status=c_prepare(int(n(2),c_int),int(n(3),c_int),int(g,c_int),v,a,s,wall)
 end subroutine
 subroutine yy_magnetic_wall_apply(ni,time,dt_old,tfac,bf,status)
  integer,intent(in)::ni
  real(c_double),intent(in)::time,dt_old,tfac
  real(c_double),contiguous,intent(inout)::bf(:,:,:,:)
  integer,intent(out)::status
  status=c_apply(wall,int(ni,c_int),time,dt_old,tfac,bf)
 end subroutine
end module
