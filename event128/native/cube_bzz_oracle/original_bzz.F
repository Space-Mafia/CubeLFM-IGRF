      SUBROUTINE BZZ(BX,BY,BZ,BI,BJ,BK)
! ... Local variables ..................................................
#include "global_dims.inc"
#include "mhd_P++_pointers.inc"
#include "fortran_P++_arrays.inc"
#include "param.inc"
#ifdef SPHERE
#include "meter.inc"
#endif
C
      COMMON /SCRACH/
     $     FDUM(NI,NJ,NK),XDUM(NIP1,NJP1,NKP1),
     $     YDUM(NIP1,NJP1,NKP1),ZDUM(NIP1,NJP1,NKP1),
     $     bxsum(nk),bysum(nk),bzsum(nk),
     $     bx_out(nk),by_out(nk),bz_out(nk)
C
      DIMENSION XAVERI(NIP1,NJ,NK),XAVERJ(NI,NJP1,NK),XAVERK(NI,NJ,NKP1)
      DIMENSION YAVERI(NIP1,NJ,NK),YAVERJ(NI,NJP1,NK),YAVERK(NI,NJ,NKP1)
      DIMENSION ZAVERI(NIP1,NJ,NK),ZAVERJ(NI,NJP1,NK),ZAVERK(NI,NJ,NKP1)
      EQUIVALENCE (XAVERI,XDUM),(XAVERJ,XDUM),(XAVERK,XDUM)
      EQUIVALENCE (YAVERI,YDUM),(YAVERJ,YDUM),(YAVERK,YDUM)
      EQUIVALENCE (ZAVERI,ZDUM),(ZAVERJ,ZDUM),(ZAVERK,ZDUM)
      DIMENSION vol_out(nk), volsum(nk)
! ... Parameter variables ..............................................
      DIMENSION BX(ILOW:IHIGH,JLOW:JHIGH,KLOW:KHIGH)
      DIMENSION BY(ILOW:IHIGH,JLOW:JHIGH,KLOW:KHIGH)
      DIMENSION BZ(ILOW:IHIGH,JLOW:JHIGH,KLOW:KHIGH)
      DIMENSION BI(ILOW:IHIGH,JLOW:JHIGH,KLOW:KHIGH)
      DIMENSION BJ(ILOW:IHIGH,JLOW:JHIGH,KLOW:KHIGH)
      DIMENSION BK(ILOW:IHIGH,JLOW:JHIGH,KLOW:KHIGH)
! ... Begin ............................................................
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in therest.F::BZZ(...)"
#endif
C
      if ( .not. MHD_PE ) return
C
      DO 200 K=1,NL_K
      DO 200 J=1,NL_J
      DO 200 I=1,NL_I
 200     FDUM(I,J,K) = BI(I,J,K)*Face_I(I,J,k)
     $        +BI(I+1,J,K)*Face_I(I+1,J,K)
C
      DO 210 K=1,NL_K
      DO 210 J=1,NL_J
      DO 210 I=1,NL_I+1
      XAVERI(I,J,K) = X(I,J,K)+X(I,J+1,K)+X(I,J,K+1)+X(I,J+1,K+1)
      YAVERI(I,J,K) = Y(I,J,K)+Y(I,J+1,K)+Y(I,J,K+1)+Y(I,J+1,K+1)
  210 ZAVERI(I,J,K) = Z(I,J,K)+Z(I,J+1,K)+Z(I,J,K+1)+Z(I,J+1,K+1)
C
      DO 220 K=1,NL_K
      DO 220 J=1,NL_J
      DO 220 I=1,NL_I
      BX(I,J,K) = FDUM(I,J,K)*(XAVERI(I+1,J,K)-XAVERI(I,J,K))
      BY(I,J,K) = FDUM(I,J,K)*(YAVERI(I+1,J,K)-YAVERI(I,J,K))
  220 BZ(I,J,K) = FDUM(I,J,K)*(ZAVERI(I+1,J,K)-ZAVERI(I,J,K))
C
C
        DO 300 K=1,NL_K
        DO 300 J=1,NL_J
        DO 300 I=1,NL_I
  300   FDUM(I,J,K)=BJ(I,J,K)*Face_J(I,J,K)+BJ(I,J+1,K)*Face_J(I,J+1,K)
C
        DO 310 K=1,NL_K
        DO 310 J=1,NL_J+1
        DO 310 I=1,NL_I
        XAVERJ(I,J,K) = X(I,J,K)+X(I+1,J,K)+X(I,J,K+1)+X(I+1,J,K+1)
        YAVERJ(I,J,K) = Y(I,J,K)+Y(I+1,J,K)+Y(I,J,K+1)+Y(I+1,J,K+1)
  310   ZAVERJ(I,J,K) = Z(I,J,K)+Z(I+1,J,K)+Z(I,J,K+1)+Z(I+1,J,K+1)
C
        DO 320 K=1,NL_K
        DO 320 J=1,NL_J
        DO 320 I=1,NL_I
        BX(I,J,K)= BX(I,J,K)+FDUM(I,J,K)*(XAVERJ(I,J+1,K)-XAVERJ(I,J,K))
        BY(I,J,K)= BY(I,J,K)+FDUM(I,J,K)*(YAVERJ(I,J+1,K)-YAVERJ(I,J,K))
  320   BZ(I,J,K)= BZ(I,J,K)+FDUM(I,J,K)*(ZAVERJ(I,J+1,K)-ZAVERJ(I,J,K))
C
C
        DO 400 K=1,NL_K
        DO 400 J=1,NL_J
        DO 400 I=1,NL_I
  400   FDUM(I,J,K)=BK(I,J,K)*Face_K(I,J,K)+BK(I,J,K+1)*Face_K(I,J,K+1)
C
        DO 410 K=1,NL_K+1
        DO 410 J=1,NL_J
        DO 410 I=1,NL_I
        XAVERK(I,J,K) = X(I,J,K)+X(I+1,J,K)+X(I,J+1,K)+X(I+1,J+1,K)
        YAVERK(I,J,K) = Y(I,J,K)+Y(I+1,J,K)+Y(I,J+1,K)+Y(I+1,J+1,K)
  410   ZAVERK(I,J,K) = Z(I,J,K)+Z(I+1,J,K)+Z(I,J+1,K)+Z(I+1,J+1,K)
C
        DO 420 K=1,NL_K
        DO 420 J=1,NL_J
        DO 420 I=1,NL_I
        BX(I,J,K)= BX(I,J,K)+FDUM(I,J,K)*(XAVERK(I,J,K+1)-XAVERK(I,J,K))
        BY(I,J,K)= BY(I,J,K)+FDUM(I,J,K)*(YAVERK(I,J,K+1)-YAVERK(I,J,K))
  420   BZ(I,J,K)= BZ(I,J,K)+FDUM(I,J,K)*(ZAVERK(I,J,K+1)-ZAVERK(I,J,K))
C
C
        DO 500 K=1,NL_K
        DO 500 J=1,NL_J
        DO 500 I=1,NL_I
        FDUM(I,J,K) = 0.125/VOLUME(I,J,K)
        BX(I,J,K) = BX(I,J,K)*FDUM(I,J,K)
        BY(I,J,K) = BY(I,J,K)*FDUM(I,J,K)
  500   BZ(I,J,K) = BZ(I,J,K)*FDUM(I,J,K)
C
#ifdef SPHERE
#ifdef DO_RINGAV
C
C    AVERAGE ALONG THE X-AXIS
C
c   need to do ringav ?
      if ( j_global_off .gt. nrings+1 .and.
     $       j_global_off + nl_j .lt. nj_global - nrings-1 ) return


c
c   loop over all j's, first for cell centered quantities
c


      do j=1,nl_j
c
c   pick front side or back side
c
         jglobal = j + j_global_off
c
#ifdef RINGAV_TEST
      write (iunit,*) 'looping over j = ',j,jglobal
#endif
c
c  front side
         if  ( jglobal .le. nrings ) then
              nring_now = jglobal
c back side
         elseif ( jglobal .ge. nj_global - nrings+1) then
            nring_now = nj_global - jglobal + 1
         else
              nring_now = 0
         endif
#ifdef RINGAV_TEST
         write (iunit, *) 'ring now ', nring_now
#endif


      if ( nring_now .gt. 0 ) then
          nav = naver(nring_now)

                
C
          do i=1,nl_i
             do k=1,nl_k
                bxsum(k) = bx(i,j,k)
                bysum(k) = by(i,j,k)
                bzsum(k) = bz(i,j,k)
             enddo
           call phi_av(bxsum,bx_out,nl_k,nring_now,nav,.true.,.false.)
           call phi_av(bysum,by_out,nl_k,nring_now,nav,.true.,.false.)
           call phi_av(bzsum,bz_out,nl_k,nring_now,nav,.true.,.false.)
             do k=1,nl_k
                bx(i,j,k) = bx_out(k)
                by(i,j,k) = by_out(k)
                bz(i,j,k) = bz_out(k)
             enddo
          enddo
       endif
       enddo
c
#endif
#endif
C

C
! ... End SUBROUTINE BZZ(BX,BY,BZ,BI,BJ,BK) .............................
        RETURN
        END
