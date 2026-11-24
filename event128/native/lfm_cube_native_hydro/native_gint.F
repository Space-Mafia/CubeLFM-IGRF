      FUNCTION GINT(FX,X,Y,Z,X1,Y1,Z1)
      SAVE
C     Integrate FX(x,y,z) over the line (X,Y,Z) to (X1,Y1,Z1).
C     This subroutine does Gaussian integration with the first twelve
C     Legendre polynomials as the basis fuctions.
C     Abromowitz and Stegun page 916.
C     
      EXTERNAL FX
      DIMENSION WT(6),A(6)
C     The A() are the positive zeros of 12th order Legendre polynomial.
      DATA A/0.1252334085,  0.3678314989,  0.5873179542,
     $       0.7699026741,  0.9041172563,  0.9815606342/
C     Next are the Gaussian Integration coefficients for a
C     12th order polynomial. 
      DATA WT/0.2491470458,  0.2334925365,  0.2031674267,
     $        0.1600783285,  0.1069393259,  0.0471753363/
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subdipol.F::GINT(...)"
#endif
      DX = (X1-X)*0.5
      DY = (Y1-Y)*0.5
      DZ = (Z1-Z)*0.5
      XBAR = 0.5*(X+X1)
      YBAR = 0.5*(Y+Y1)
      ZBAR = 0.5*(Z+Z1)
C
      SUM = 0.
C
      DO 200 K=1,6
 200     SUM = SUM + WT(K)*(
     $        FX(XBAR+A(K)*DX,YBAR+A(K)*DY,ZBAR+A(K)*DZ) +
     $        FX(XBAR-A(K)*DX,YBAR-A(K)*DY,ZBAR-A(K)*DZ) )
C
      GINT = 0.5*SUM
C
      RETURN
      END
