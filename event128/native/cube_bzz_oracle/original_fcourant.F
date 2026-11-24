        subroutine fcourant(crnt_local)
! ... Variables ........................................................
#include "global_dims.inc"
#include "mhd_P++_pointers.inc"
#include "fortran_P++_arrays.inc"
#include "param.inc"
#include "help.inc"
#include "bzero.inc"
***#include "var.inc"
***#include "var1.inc"
#include "var2.inc"
***#include "meter.inc"
#include "twod.inc"
#include "run-time.inc"
c
c
        DIMENSION BSQ(NI,NJ,NK),DUM(NI,NJ,NK),V(NI,NJ,NK),
     $          VSIG(NI,NJ,NK),dtq(ni,nj,nk),
     $          dtqi(ni,nj,nk),dtqj(ni,nj,nk),dtqk(ni,nj,nk)
        EQUIVALENCE (BSQ,RHO2),(DUM,VX2),(V,VY2),(VSIG,VZ2),(dtq,c2),
     $          (dtqi,bi2),(dtqj,bj2),(dtqk,bk2)
        dimension cycle(ni,4)
! Read-only diagnostic payload; no change to the CFL expression or ordering.
        real cfl_diag(18)
        common /lfm_cfl_diagnostic/ cfl_diag
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in fcourant.F::fcourant(crnt_local)"
#endif
C
        if ( .not. MHD_PE ) return
C
        DO 210 i=1,Nl_i
        DO 210 K=1,Nl_k
        DO 210 J=1,Nl_J
        BSQ(i,j,k) = (BX_g(i,j,k)**2 + BY_g(i,j,k)**2 + BZ_g(I,J,K)**2
     $         +0.5*(bsqqk(i,j,k)+bsqqk(i,j,k+1)) )*PI4INV
        DUM(i,j,k) = 1./( RHO_g(i,j,k) + BSQ(i,j,k)*CAINV**2 )
  210   CONTINUE
C
C
        DO 300 K=1,NL_K
        DO 300 J=1,NL_J
        DO 300 I=1,NL_I
        V(I,J,K) = SQRT(VX_g(I,J,K)**2+VY_g(I,J,K)**2+VZ_g(I,J,K)**2)
        VSIG(i,j,k) = SQRT( C_g(i,j,k)**2 +
     $           BSQ(i,j,k)*DUM(i,j,k)) + V(i,j,k)
  300 CONTINUE
C
      crnt_local = 1.e30
      do k=1,nl_k
         do j=1,nl_j
            do i=1,nl_i
#ifndef BUGCOURANT
               crnt_local = min(crnt_local,
     $              cell_length(i,j,k)/vsig(i,j,k))
               if (cell_length(i,j,k)/vsig(i,j,k)
     $              .eq. crnt_local) then
                  ii0=i
                  jj0=j
                  kk0=k
               endif
#else
               ccc_0 = cell_length(I,j,k)/vsig(i,j,k)
               if ( ccc_0 .lt. crnt_local ) then
                  ii0 = i
                  jj0 = j
                  kk0 = k
                  crnt_local = ccc_0
               endif
#endif
            enddo
         enddo
      enddo
      cfl_diag(1)=crnt_local
      cfl_diag(2)=ii0+I_GLOBAL_OFF
      cfl_diag(3)=jj0+J_GLOBAL_OFF
      cfl_diag(4)=kk0+K_GLOBAL_OFF
      cfl_diag(5)=cell_length(ii0,jj0,kk0)
      cfl_diag(6)=vsig(ii0,jj0,kk0)
      cfl_diag(7)=rho_g(ii0,jj0,kk0)
      cfl_diag(8)=c_g(ii0,jj0,kk0)
      cfl_diag(9)=vx_g(ii0,jj0,kk0)
      cfl_diag(10)=vy_g(ii0,jj0,kk0)
      cfl_diag(11)=vz_g(ii0,jj0,kk0)
      cfl_diag(12)=bx_g(ii0,jj0,kk0)
      cfl_diag(13)=by_g(ii0,jj0,kk0)
      cfl_diag(14)=bz_g(ii0,jj0,kk0)
      cfl_diag(15)=bsq(ii0,jj0,kk0)
      cfl_diag(16)=bsqqk(ii0,jj0,kk0)
      cfl_diag(17)=v(ii0,jj0,kk0)
      cfl_diag(18)=dum(ii0,jj0,kk0)
C
#ifdef BUGCOURANT
      write (6,*) 'At processor #: ', mype, ii0, jj0, kk0, crnt_local
      write (6,911)
     $     mype, VSIG(ii0,jj0,kk0), V(ii0,jj0,kk0),
     $     Bsq(ii0,jj0,kk0),dum(ii0,jj0,kk0)
      write (6,912)
     $     mype, rho_g(ii0,jj0,kk0),vx_g(ii0,jj0,kk0),
     $     vy_g(ii0,jj0,kk0),vz_g(ii0,jj0,kk0)
      write (6,913)
     $     mype, c_g(ii0,jj0,kk0),bx_g(ii0,jj0,kk0),
     $     by_g(ii0,jj0,kk0),bz_g(ii0,jj0,kk0)
 911  format(i4,' -- Vs,V,Bsq,dum: ',2x,4(1pe12.4,2x))
 912  format(i4,' -- rho,vxyz: ',2x,4(1pe12.4,2x))
 913  format(i4,' -- c, bxyz:  ',2x,4(1pe12.4,2x))
#endif

!... End subroutine fcourant(crnt_local)...............................
      return
      end
