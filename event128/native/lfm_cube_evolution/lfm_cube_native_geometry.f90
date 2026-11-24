! Array adapter only: actual original METRIC/BZZ/FCOURANT/FLUX_FIX are linked.
module lfm_cube_native_geometry
 use iso_fortran_env,only:real64
 use lfm_cube_mesh
 use lfm_cube_native_storage
 implicit none
 private
 public::native_load,native_recover,native_cfl,native_ct
contains
 subroutine native_load(p)
  type(cube_patch),intent(in)::p
  integer::n(3)
  n=p%layout%local
  if(any(n>[64,32,32])) error stop 'native local capacity exceeded'
  nl_i=n(1);nl_j=n(2);nl_k=n(3)
  x=0;y=0;z=0
  x(1:n(1)+1,1:n(2)+1,1:n(3)+1)=p%x(1,:,:,:)
  y(1:n(1)+1,1:n(2)+1,1:n(3)+1)=p%x(2,:,:,:)
  z(1:n(1)+1,1:n(2)+1,1:n(3)+1)=p%x(3,:,:,:)
  call lfm_legacy_scratch_ready()
  call metric()
  cell_length(1:n(1),1:n(2),1:n(3))=p%cfl_length
 end subroutine
 subroutine native_recover(p,flux,correction,b)
  type(cube_patch),intent(in)::p
  real(real64),intent(in)::flux(0:,0:,0:,:),correction(:,:,:,:,:)
  real(real64),intent(out)::b(:,:,:,:)
  integer::i,j,k,n(3)
  ! BZZ reads only the staggered active slabs overwritten below.
  ! Do not clear unused fixed-capacity padding on every small MPI tile.
  n=p%layout%local
  bi_g(1:n(1)+1,1:n(2),1:n(3))=flux(0:n(1),1:n(2),1:n(3),1)/face_i(1:n(1)+1,1:n(2),1:n(3))
  bj_g(1:n(1),1:n(2)+1,1:n(3))=flux(1:n(1),0:n(2),1:n(3),2)/face_j(1:n(1),1:n(2)+1,1:n(3))
  bk_g(1:n(1),1:n(2),1:n(3)+1)=flux(1:n(1),1:n(2),0:n(3),3)/face_k(1:n(1),1:n(2),1:n(3)+1)
  call bzz(bx_g,by_g,bz_g,bi_g,bj_g,bk_g)
  do k=1,n(3)
   do j=1,n(2)
    do i=1,n(1)
     b(i,j,k,:)=matmul(correction(:,:,i,j,k),[bx_g(i,j,k),by_g(i,j,k),bz_g(i,j,k)])
    enddo
   enddo
  enddo
 end subroutine
 function native_cfl(q,n,ca,background_sq) result(dt)
  integer,intent(in)::n(3)
  real(real64),intent(in)::q(:,:,:,:),ca
  real(real64),intent(in)::background_sq(0:,0:,0:)
  real(real64)::dt
  rho_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,1);c_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,5)
  vx_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,2);vy_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,3)
  vz_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,4)
  bx_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,6);by_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,7)
  bz_g(1:n(1),1:n(2),1:n(3))=q(:,:,:,8)
  bsqqk(1:n(1),1:n(2),1:n(3)+1)=background_sq(1:n(1),1:n(2),0:n(3))
  pi4inv=1/(4*acos(-1._real64));cainv=1/ca
  call fcourant(dt)
 end function
 subroutine native_ct(flux,emf,n,dt)
  integer,intent(in)::n(3)
  real(real64),intent(inout)::flux(0:,0:,0:,:)
  real(real64),intent(in)::emf(0:,0:,0:,:),dt
  ! Original native E has the v-cross-B sign, unlike physical E=-v-cross-B.
  ! Exact FLUX_FIX/AVE2FLUX/FLUX2AVE use only the active face/edge
  ! slabs below, each fully overwritten. Unread capacity is left alone.
  bi_dbl(1:n(1)+1,1:n(2),1:n(3))=flux(:,1:n(2),1:n(3),1)
  bj_dbl(1:n(1),1:n(2)+1,1:n(3))=flux(1:n(1),:,1:n(3),2)
  bk_dbl(1:n(1),1:n(2),1:n(3)+1)=flux(1:n(1),1:n(2),:,3)
  ei_g(1:n(1),1:n(2)+1,1:n(3)+1)=emf(1:n(1),:,:,1)
  ej_g(1:n(1)+1,1:n(2),1:n(3)+1)=emf(:,1:n(2),:,2)
  ek_g(1:n(1)+1,1:n(2)+1,1:n(3))=emf(:,:,1:n(3),3)
  call flux_fix(dt)
  flux(:,1:n(2),1:n(3),1)=bi_dbl(1:n(1)+1,1:n(2),1:n(3))
  flux(1:n(1),:,1:n(3),2)=bj_dbl(1:n(1),1:n(2)+1,1:n(3))
  flux(1:n(1),1:n(2),:,3)=bk_dbl(1:n(1),1:n(2),1:n(3)+1)
 end subroutine
end module
