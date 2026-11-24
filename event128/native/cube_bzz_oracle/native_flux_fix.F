      subroutine flux_fix(DT)
! ... Local Variables ..................................................
#include "global_dims.inc"
#include "mhd_P++_pointers.inc"
#include "fortran_P++_arrays.inc"
#include "param.inc"
#include "help.inc"
#include "var.inc"
#ifdef BZEROS
#include "bzero.inc"
#endif
#include "b_dbl.inc"
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in hack.F::flux_fix(DT)"
#endif
      nl_ip1 = nl_i+1
      nl_jp1 = nl_j+1
      nl_kp1 = nl_k+1
c
      do k=1,nl_k
         do j=1,nl_j
            do i=1,nl_ip1
               bi_dbl(i,j,k) = bi_dbl(i,j,k) + DBLE(dt)*(
     $               DBLE(Ek_g(i,j+1,k)) - DBLE(Ek_g(i,j,k))
     $              -DBLE(Ej_g(i,j,k+1))  + DBLE(Ej_g(i,j,k)) )
            enddo
         enddo
      enddo
            
c
      do k=1,nl_k
         do j=1,nl_jp1
            do i=1,nl_i
               bj_dbl(i,j,k) = bj_dbl(i,j,k) + DBLE(dt)*(
     $               DBLE(Ei_g(i,j,k+1)) - DBLE(Ei_g(i,j,k))
     $              -DBLE(Ek_g(i+1,j,k))  + DBLE(Ek_g(i,j,k)) )
            enddo
         enddo
      enddo
            
c
      do k=1,nl_kp1
         do j=1,nl_j
            do i=1,nl_i
               bk_dbl(i,j,k) = bk_dbl(i,j,k) + DBLE(dt)*(
     $               DBLE(Ej_g(i+1,j,k)) - DBLE(Ej_g(i,j,k))
     $              -DBLE(Ei_g(i,j+1,k))  + DBLE(Ei_g(i,j,k)) )
            enddo
         enddo
      enddo
c
c
      call ave2flux
c
      do k=1,nl_kp1
         do j = 1,nl_j
            do i=1,nl_i
               bk_g(i,j,k) = bk_dbl(i,j,k)
            enddo
         enddo
      enddo
      do k=1,nl_k
         do j = 1,nl_jp1
            do i=1,nl_i
               bj_g(i,j,k) = bj_dbl(i,j,k)
            enddo
         enddo
      enddo
      do k=1,nl_k
         do j = 1,nl_j
            do i=1,nl_ip1
               bi_g(i,j,k) = bi_dbl(i,j,k)
            enddo
         enddo
      enddo
c
      call flux2ave
c
#ifdef BZEROS
      do k=1,nl_kp1
         do j = 1,nl_j
            do i=1,nl_i
               bk_g(i,j,k) = bk_g(i,j,k)-bnqfk(i,j,k)
            enddo
         enddo
      enddo
      do k=1,nl_k
         do j = 1,nl_jp1
            do i=1,nl_i
               bj_g(i,j,k) = bj_g(i,j,k)-bnqfj(i,j,k)
            enddo
         enddo
      enddo
      do k=1,nl_k
         do j = 1,nl_j
            do i=1,nl_ip1
               bi_g(i,j,k) = bi_g(i,j,k)-bnqfi(i,j,k)
            enddo
         enddo
      enddo
#endif
!
! ... End subroutine FLUX_FIX(DT) ......................................
!
      return
      end
