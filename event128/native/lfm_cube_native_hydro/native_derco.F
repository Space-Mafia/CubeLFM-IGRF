      SUBROUTINE DERCO(FCOEF,NORDER, local_proc)
! ... Local variables ..................................................
c      IF ( LOCAL_PROC .EQ. 0) WRITE (9,112) ORDER,(FCOEF(I),I=1,ORDER)
c      NORDER-TH ORDER CALCULATION for the fluxes F use
c      CENTERED DERIVATIVES USE THE FOLLOWING COEFFICIENTS
c
      DIMENSION FCOEF(NORDER)
      COMMON /SCRACH/
     $     A(200),B(200),C(200),F(200)
! ... Begin ............................................................
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in init-fortran.F::DERCO(...)"
#endif
C
      N = NORDER
      N2 = N/2
C
      A(1) = 1.
      DO 200 I = 2,N2
  200 A(I) = A(I-1)*FLOAT(N+2-I)/FLOAT(I-1)
C
      DO 201 I = 1,N2
  201 B(I) = A(N2+1-I)
C
      DO 202 I = 1,N2,2
  202 B(I) = B(I)/FLOAT(I)
C
      DO 203 I = 2,N2,2
  203 B(I) = - B(I)/FLOAT(I)
C
      IF ( LOCAL_PROC .EQ. 0) WRITE (9,*)(A(I),I=1,N2)
  101 FORMAT (5(1PE16.4))
      IF ( LOCAL_PROC .EQ. 0) WRITE (9,*) (B(I),I=1,N2)
      C(1) = 0.
      DO 300 I = 1,N2
  300 C(I+1) = B(I)
C
      F(1) = B(N2)
      DO 400 I = 1,N2
  400 F(I+1) = C(N2+1-I) + F(I)
C
      N2P1 = N2+1
      IF ( LOCAL_PROC .EQ. 0) WRITE(9,*) (C(I),I=1,N2P1)
      IF ( LOCAL_PROC .EQ. 0) WRITE (9,*) (F(I),I=1,N2P1)
      SUM = 0.
      DO 401 I = 1,N2
  401 SUM = SUM + F(I)
      SUM = SUM*2.
C
      DO 410 I = 1,N2
  410 F(I) = F(I)/SUM
C
      DO 500 I = 1,N2
      FCOEF(I) = F(I)
  500 FCOEF(N+1-I) = F(I)
!
! ... End SUBROUTINE DERCO ..............................................
!
      RETURN
      END
