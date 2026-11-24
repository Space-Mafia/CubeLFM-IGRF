C**********************************************************************
      SUBROUTINE Z2ICG(FI,A1,A2,A3,A4,A5,A6,A7,A8,A9,SI,NS5,
     .                 NX,NY,NDX,NDY,LBOUN,ITMAX,ITER,ERROR,SER,
     .                  AB1,s1,AB2,s2)
      SAVE
C
C This module solves a general system of equations of the form
C
C        + A1(J,K)*FI(J-1,K  )
C        + A2(J,K)*FI(J  ,K-1)
C        + A3(J,K)*FI(J,K)
C        + A4(J,K)*FI(J  ,K+1)
C        + A5(J,K)*FI(J+1,K  )
C        + A6(J,K)*FI(J-1,K-1) + A7(J,K)*FI(J-1,K+1)
C        + A8(J,K)*FI(J+1,K-1) + A9(J,K)*FI(J+1,K+1) = SI(J,K)
C
C for J = 2, NX-1 and K = 2, NY-1.
C
C The stencil surrounding a point (J,K) looks like
C
C                7 | 4 | 9        K ^
C                ----------         |
C                1 | 3 | 5          |
C                ----------         --------->
C                6 | 2 | 8                   J
C
C The cells (J=1;K=1,NY), (J=NX;K=1,NY), (J=1,NX;K=1), (J=1,NX;K=NY)
C are boundary cells.
C
C An outline of the code is as follows:
C
C  (i) The system of eqns are manipulated for B.C.:
C      For Dirichlet on a side the coefficients A multiplying
C      boundary terms of F are redefined so that they are
C      included in the source term.
C      For Neumann B.C. on a side the coefficients A multiplying
C      boundary terms of F are incorporated into other coefficients.
C      For periodic B.C. on a side the coefficients A multiplying
C      boundary terms of F are incorporated into the main matrix
C      as the calculation proceeds.
C
C (ii) The five point stencil linking points 1,2,3,4,5 is split into a
C      symmetric system and a asymmetric one according to
C
C        AS1(J,K)*F(J-1,K)   + AA1(J,K)*F(J-1,K)
C      + AS2(J,K)*F(J,K-1)   + AA2(J,K)*F(J,K-1)
C      + AS3(J,K)*F(J,K)
C      + AS2(J,K+1)*F(J,K+1) - AS2(J,K+1)*F(J,K+1)
C      + AS1(J+1,K)*F(J+1,K) - AS1(J+1,K)*F(J+1,K)
C      + A6(J,K)*F(J-1,K-1)  + A7(J,K)*F(J-1,K+1)
C      + A8(J,K)*F(J+1,K-1)  + A9(J,K)*F(J+1,K+1) = S(J,K)
C
C (iii) The symmetric part is solved iteratively using an ICCG solver
C       through Z2SOL, Z2SYM, Z2INV, and Z2MAT.
C       (Incomplete Cholesky (or Crout) decomposition with
C        Conjugate Gradient technique.)
C       Some references are:
C          i) D. S. Kershaw, J. Comp. Phys., 26, 43, (1978),
C             "The Incomplete Cholesky-Conjugate Gradient method for
C              the Iterative solution of Systems of Equations."
C         ii) G. H. Golub and C. F. Van Lear, Matrix Computaions,
C             (Johns Hopkins Univerisity Press: Baltimore),
C              1983, chapter 10.
C        iii) P. Concus, G. H. Golub, and D. P. O'Leary, "A
C             Generalized Conjugate Gradient Method for the Numerical
C             Solution of Elliptic Partial Differential Equations,"
C             in Sparse Matrix Computations, ed. J. R. Bunch and
C             D. J. Rose (Academic press: New York), 1976, p.304.
C             [simplest proof of ICCG technique).
C         iv) D. M. Young and K. C. Jea, "Generalized Conjugate
C             Gradient Acceleration on Non-symmetrizable Iterative
C             Methods," in Large scale Matrix Problems, ed. Ake
C             Bjorck, R. J. Plemmons, and H. Schneider, (North
C             Holland: New York), 1981, p.159.
C             [possible techniques for non-symmetric terms]
C
C (iv) Given the solution of the symmetric 5-point system the
C      asymmetric and other terms are solved via an iteration
C      scheme in Z2ASY.  The scheme is based on the iteration formula
C
C               l+1    l          l    l   l-1
C              F    = F  + w ( a Z  + F - F    )
C                           l   l
C
C      and the residual defines Z
C
C                    l          l
C                 M*Z  = S - A*F
C                                                    T
C      where l is the iteration number and M = (L*D*L ), i.e., the
C      incomplete decomposition of the 5 point symmetric system.
C      These formulas are from Ref(iii). In this ref the
C      parameters a and w are determined so that
C             l     l+1       l-1     l+1
C           <Z , M*Z   > = < Z   , M*Z   > = 0
C      which leads to the conjugate gradient scheme as long as A
C      is symmetric. In the present scheme, where A may be
C      non-symmetric, the parameter a is determined such that
C             l     l+1
C           <Z , M*Z   > = 0
C      and w is determined such that the norm of the residual
C               l+1     l+1
C           <M*Z   , M*Z   >
C       is minimized.  It has been found in practice that this
C       iteration scheme does not insure convergence to an
C       arbitrarily small error, but of all the methods
C       investigated this one reduces the residual norm significantly
C       in the fewest iterations.  Other approaches, such as
C       Tchebychev iteration for non-symmetric systems, require
C       a good knowledge of A's eigenvalues and reqiure all
C       eigenvalues to have no zero real parts. These conditions
C       are uncertain for arbitrary systems intended to be used here.
C
C The general flow of the code is modeled after Klaus Hain's ICCG
C potential solver except that recursive procedures are used to
C determine the LDU decomposition and invert triangular matrices.
C There are changes in the form of the matrix eqn depending
C on the boundary conditions.  This is done to speed
C convergence for Neuman Boundary conditions.
C There are also changes in the iterative solver for
C asymmetric terms and the additional four terms of a nine point
C scheme.
C
C       FI    = Initial value of unknown vector on input,
C             = on output it is the new solution (OUTPUT).
C       A'S   = Coefficients of the equation to be solved.
C               (note: the A's might be changed in this module.)
C       SI    = Initial source term.
C               (note: the SI's might be changed in this module.)
C       NS5   = 1, if one has a purely 5 point symmetric system,
C               0,  otherwise.
C       NX    = Range in J index to be solved
C               including two boundary points.
C       NY    = Range in K index to be solved
C               including two boundary points.
C       NDX   = J dimension of incoming arrays.
C       NDY   = K dimension of incoming arrays.
C               (note: presently NDX<52, NDY<52, but can be
C                changed through MX and MY in PARAMETER statements
C                in each module).
C       LBOUN(I) = Specifies boundary conditions on each side.
C                                                     (I=4)  K=NY
C                  I = 1 B. C. at J = 1  (XA)      |--------|
C                  I = 2 B. C. at J = NX (XE) (I=1)|        | (I=2)
C                  I = 3 B. C. at K = 1  (YA)      |        |
C                  I = 4 B. C. at K = NY (YE)      |--------| K=1
C                                               J=1  (I=3)  J=NX
C                  For e.g.,
C                  LBOUN(1) = 0 Dirichlet          : F(1,K) specified.
C                  LBOUN(1) = 1 Neumann(DF/DN = 0) : F(1,K) = F(2,K)
C                  LBOUN(1) = 2 anti-Neuman        : F(1,K) =-F(2,K)
C                  LBOUN(1) = 3 periodic           : F(1,K) = F(NX-1,K)
C                  LBOUN(1) = 4 odd periodic       : F(1,K) =-F(NX-1,K)
C                and similarly for I = 2,3,4.
C
C       ITMAX    = Maximum number of iterations allowed.
C       ITER     = Number of iterations used (OUTPUT).
C       ERROR    = Input error condition for convergence.
C       SER      = Error of output solution (OUTPUT).
C DATA statement within code contains
C       LPRINT = if .TRUE. prints out convergence errors AT
C       each iteration.
C
C       This module calls the subroutines and entries
C       ZSOLV, Z2SYM, Z2ASY, Z2INV, Z2MAT.
C
C----------------------------------------------------------------------
C......dimension block
C
#include "z2param.inc"
C
      REAL FI(NDX,NDY),
     .     A1(NDX,NDY),A2(NDX,NDY),A3(NDX,NDY),A4(NDX,NDY),A5(NDX,NDY),
     .     A6(NDX,NDY),A7(NDX,NDY),A8(NDX,NDY),A9(NDX,NDY),SI(NDX,NDY),
     .     AB1(NDY),AB2(NDY)
C
*      DOUBLE PRECISION F,AS1,AA1,AS2,AA2,AS3,S,SA,PZ,R,P,Q,H,Z,G    !VAX
C
      COMMON/ICCG1/ AS1(MX,MY),AA1(MX,MY),AS2(MX,MY),AA2(MX,MY),
     .              AS3(MX,MY),PZ(MX,MY),H(MX,MY),Z(MX,MY),G(MX,MY),
     .              NX0,NY0,NX1,NY1,LPRINT
      COMMON/ICCG2/ R(MX,MY),P(MX,MY),Q(MX,MY),
     .              F(MX,MY),S(MX,MY),SA(MX,MY)
C
      INTEGER LBOUN(4)
      LOGICAL LPRINT
#ifdef DEBUG
      DATA    LPRINT/.true./
#else
      DATA    LPRINT/.false./
#endif 
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in z2code.F::Z2ICG(...)"
#endif
C
C----------------------------------------------------------------------
      NX0 = NX
      NX1 = NX-1
      NY0 = NY
      NY1 = NY-1
C----------------------------------------------------------------------
C.....redefine A's for Dirichlet B.C.
      IF(LBOUN(1) .EQ. 0) THEN
        DO 1 K = 2, NY1
          SI(2,K) = SI(2,K)-A1(2,K)*FI(1,K)
          A1(2,K) = 0.0
1       CONTINUE
      ENDIF
      IF(LBOUN(2) .EQ. 0) THEN
        DO 2 K = 2, NY1
          SI(NX1,K) = SI(NX1,K)-A5(NX1,K)*FI(NX,K)
          A5(NX1,K) = 0.0
2       CONTINUE
      ENDIF
      IF(LBOUN(3) .EQ. 0) THEN
        DO 3 J = 2, NX1
          SI(J,2) = SI(J,2)-A2(J,2)*FI(J,1)
          A2(J,2) = 0.0
*  try to fix for finite element case
          a8(j,2) = 0.0
          a6(j,2) = 0.0
3       CONTINUE
      ENDIF
      IF(LBOUN(4) .EQ. 0) THEN
        DO 4 J = 2, NX1
          SI(J,NY1) = SI(J,NY1)-A4(J,NY1)*FI(J,NY)
          A4(J,NY1) = 0.0
*  try to fix for finite element case
          a7(j,ny1) = 0.0
          a9(j,ny1) = 0.0
4       CONTINUE
      ENDIF
C
C.....rearrange A's for Neumann B.C.:
      IF(LBOUN(1) .EQ. 1 .OR. LBOUN(1) .EQ. 2) THEN
        SGN = 3-2*LBOUN(1)
        DO 5 K = 2, NY1
          A3(2,K) = A3(2,K)+SGN*A1(2,K)
          A1(2,K) = 0.0
5       CONTINUE
      ENDIF
      IF(LBOUN(2) .EQ. 1 .OR. LBOUN(2) .EQ. 2) THEN
        SGN = 3-2*LBOUN(2)
        DO 6 K = 2, NY1
          A3(NX1,K) = A3(NX1,K)+SGN*A5(NX1,K)
          A5(NX1,K) = 0.0
6       CONTINUE
      ENDIF
      IF(LBOUN(3) .EQ. 1 .OR. LBOUN(3) .EQ. 2) THEN
        SGN = 3-2*LBOUN(3)
        DO 7 J = 2, NX1
          A3(J,2) = A3(J,2)+SGN*A2(J,2)
          A2(J,2) = 0.0
7       CONTINUE
      ENDIF
      IF(LBOUN(4) .EQ. 1 .OR. LBOUN(4) .EQ. 2) THEN
        SGN = 3-2*LBOUN(4)
        DO 8 J = 2, NX1
          A3(J,NY1) = A3(J,NY1)+SGN*A4(J,NY1)
          A4(J,NY1) = 0.0
8       CONTINUE
      ENDIF
C
C----------------------------------------------------------------------
C.....form symmetric and antisymmetric 5 point subsysytem:
      DO 21 K = 2, NY1
      DO 21 J = 2, NX1
        AS3(J,K) = A3(J,K)
21    CONTINUE
      DO 23 K = 2, NY1
        IF(LBOUN(1) .GE. 3) THEN
          AS1(2,K)  = 0.5*(A1(2,K)+A5(NX1,K))
          AA1(2,K)  = 0.5*(A1(2,K)-A5(NX1,K))
        ELSE
          AS1(2,K)  = 0.0
          AA1(2,K)  = 0.0
        ENDIF
        AS1(NX,K) = AS1(2,K)
        AA1(NX,K) =-AA1(2,K)
        DO 22 J = 3, NX1
          AS1(J,K) = 0.5*(A1(J,K)+A5(J-1,K))
          AA1(J,K) = 0.5*(A1(J,K)-A5(J-1,K))
22      CONTINUE
23    CONTINUE
      DO 25 J = 2, NX1
        IF(LBOUN(3) .GE. 3) THEN
          AS2(J,2)  = 0.5*(A2(J,2)+A4(J,NY1))
          AA2(J,2)  = 0.5*(A2(J,2)-A4(J,NY1))
        ELSE
          AS2(J,2)  = 0.0
          AA2(J,2)  = 0.0
        ENDIF
        AS2(J,NY) = AS2(J,2)
        AA2(J,NY) =-AA2(J,2)
        DO 24 K = 3, NY1
          AS2(J,K) = 0.5*(A2(J,K)+A4(J,K-1))
          AA2(J,K) = 0.5*(A2(J,K)-A4(J,K-1))
24      CONTINUE
25    CONTINUE
C
C----------------------------------------------------------------------
C.....transfer to double precision:
      DO 31 K = 1, NY
      DO 31 J = 1, NX
        F(J,K)  = FI(J,K)
31    CONTINUE
      DO 32 K = 2, NY1
      DO 32 J = 2, NX1
        S(J,K)  = SI(J,K)
32    CONTINUE
C
C.....now solve matrix problem:
      CALL Z2SOL(NS5,LBOUN,ITMAX,ITER,ERROR,SER,A6,A7,A8,A9,NDX,NDY,
     .                  AB1,s1,AB2,s2)
C
C.....put back solution in single precision:
      DO 33 K = 1, NY
      DO 33 J = 1, NX
        FI(J,K) = F(J,K)
33    CONTINUE
C
      RETURN
      END
C**********************************************************************
C**********************************************************************
C**********************************************************************
      SUBROUTINE Z2INV (Q,R,LBOUN)
      SAVE
C
C.....dimension block:
C
#include "z2param.inc"
C
*      DOUBLE PRECISION AS1,AA1,AS2,AA2,AS3,S,PZ,R,Q,H,Z,G    !VAX
C
      COMMON/ICCG1/ AS1(MX,MY),AA1(MX,MY),AS2(MX,MY),AA2(MX,MY),
     .              AS3(MX,MY),PZ(MX,MY),H(MX,MY),Z(MX,MY),G(MX,MY),
     .              NX,NY,NX1,NY1,LPRINT
      DIMENSION R(MX,MY),Q(MX,MY),S(MX,MY)
C
      INTEGER LBOUN(4)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in z2code.F::Z2INV(...)"
#endif
C----------------------------------------------------------------------
C.....do an L*D*LT inversion:
C
      DO 1 K = 2, NY1
        H(1,K)  = 0.0
        Q(NX,K) = 0.0
1     CONTINUE
      DO 2 J = 2, NX1
        H(J,1)  = 0.
        Q(J,NY)  = 0.
2     CONTINUE
C
      FX = (LBOUN(1)/3)*(1 - 2*(LBOUN(1)/4))
      FY = (LBOUN(3)/3)*(1 - 2*(LBOUN(3)/4))
C
C.....do forward solution:
      DO 12 J = 2, NX1-1
        DO 11 K = 2, NY1
          H(J,K) = PZ(J,K)*R(J,K)-PZ(J,K)*AS1(J,K)*H(J-1,K)
     .                           -PZ(J,K)*AS2(J,K)*H(J,K-1)
11      CONTINUE
        H(J,NY1) = H(J,NY1)-FY*PZ(J,NY1)*AS2(J,2)*H(J,2)
12    CONTINUE
      DO 13 K = 2, NY1
        H(NX1,K) = PZ(NX1,K)*R(NX1,K)-PZ(NX1,K)*AS1(NX1,K)*H(NX1-1,K)
     .                               -PZ(NX1,K)*AS2(NX1,K)*H(NX1,K-1)
     .                               -FX*PZ(NX1,K)*AS1(2,K)*H(2,K)
13    CONTINUE
      H(NX1,NY1) = H(NX1,NY1)-FY*PZ(NX1,NY1)*AS2(NX1,2)*H(NX1,2)
C
C.....now backward substitution:
      DO 15 J = NX1, 3, -1
        DO 14 K = NY1, 2, -1
          Q(J,K) = H(J,K)-PZ(J,K)*AS2(J,K+1)*Q(J,K+1)
     .                   -PZ(J,K)*AS1(J+1,K)*Q(J+1,K)
14      CONTINUE
        Q(J,2) = Q(J,2)-FY*PZ(J,2)*AS2(J,2)*Q(J,NY1)
15    CONTINUE
      DO 16 K = NY1, 2, -1
        Q(2,K) = H(2,K)-PZ(2,K)*AS2(2,K+1)*Q(2,K+1)
     .                 -PZ(2,K)*AS1(3,K)*Q(3,K)
     .                 -FX*PZ(2,K)*AS1(2,K)*Q(NX1,K)
16    CONTINUE
      Q(2,2) = Q(2,2)-FY*PZ(2,2)*AS2(2,2)*Q(2,NY1)
C
      IF ( LBOUN(1) .GE .3 .AND. LBOUN(3) .GE. 3) THEN
        DO 21 K=2,NY1
        DO 21 J=2,NX1
          Q(J,K) = Q(J,K) - Q(2,2)
21      CONTINUE
      ENDIF
C
      RETURN
C
C**********************************************************************
      ENTRY Z2MAT (Q,R,S,NS,LBOUN)
C
      IF ( LBOUN(1).GE.3 ) THEN
        FX = 1 - 2*(LBOUN(1)/4)
        DO 32 K=2,NY1
          R(1,K)  = FX*R(NX1,K)
          R(NX,K) = FX*R(2,K)
32      CONTINUE
      ENDIF
      IF ( LBOUN(3).GE.3 ) THEN
        FY = 1 - 2*(LBOUN(3)/4)
        DO 33 J=2,NX1
          R(J,1)  = FY*R(J,NY1)
          R(J,NY) = FY*R(J,2)
33      CONTINUE
      ENDIF
C
      DO 34 K=2,NY1
      DO 34 J=2,NX1
        Q(J,K) =   AS1(J,K)*R(J-1,K)+AS2(J,K)*R(J,K-1)
     .            +AS3(J,K)*R(J,K)
     .            +AS2(J,K+1)*R(J,K+1)+AS1(J+1,K)*R(J+1,K)
34    CONTINUE
C
      IF(NS .EQ. 2) RETURN
C
      DO 35 K=2,NY1
      DO 35 J=2,NX1
        Q(J,K) = S(J,K) - Q(J,K)
35    CONTINUE
C
      RETURN
C
      END
C**********************************************************************
C**********************************************************************
C**********************************************************************
      SUBROUTINE Z2SOL(NS5,LBOUN,ITMAX,ITER,ERROR,SER,
     .           A6,A7,A8,A9,NDX,NDY,AB1,s1,AB2,s2)
      SAVE
C
C......dimension block:
C
      REAL A6(NDX,NDY),A7(NDX,NDY),A8(NDX,NDY),A9(NDX,NDY)
        REAL AB1(NDY),AB2(NDY)
C
#include "z2param.inc"
C
*      DOUBLE PRECISION F,AS1,AA1,AS2,AA2,AS3,S,SA,PZ,R,P,Q,H,Z,G,    !VAX
*     .                 RS,QS                                        !VAX
C
      COMMON/ICCG1/ AS1(MX,MY),AA1(MX,MY),AS2(MX,MY),AA2(MX,MY),
     .              AS3(MX,MY),PZ(MX,MY),H(MX,MY),Z(MX,MY),G(MX,MY),
     .              NX,NY,NX1,NY1,LPRINT
      COMMON/ICCG2/ R(MX,MY),P(MX,MY),Q(MX,MY),
     .              F(MX,MY),S(MX,MY),SA(MX,MY)
C
      INTEGER LBOUN(4)
      LOGICAL LPRINT
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in z2code.F::Z2SOL(...)"
#endif
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C......find diagonal elements of L*D*LT decomposition by recursion:
C
      DO 1 K=2,NY1
        PZ(1,K) = 0.0
1     CONTINUE
      DO 2 J = 2, NX1
        PZ(J,1) = 0.0
2     CONTINUE
C
      FX = LBOUN(1)/3
      FY = LBOUN(3)/3
C
      DO 12 J = 2, NX1-1
        DO 11 K = 2, NY1-1
          PZ(J,K) = 1./(AS3(J,K)-AS1(J,K)*PZ(J-1,K)*AS1(J,K)
     .                          -AS2(J,K)*PZ(J,K-1)*AS2(J,K) )
11      CONTINUE
        PZ(J,NY1) = 1./(AS3(J,NY1)-AS1(J,NY1)*PZ(J-1,NY1)*AS1(J,NY1)
     .                            -AS2(J,NY1)*PZ(J,NY1-1)*AS2(J,NY1)
     .                            -FY*AS2(J,2)*PZ(J,2)*AS2(J,2)      )
12    CONTINUE
      DO 13 K = 2, NY1-1
        PZ(NX1,K) = 1./(AS3(NX1,K)-AS1(NX1,K)*PZ(NX1-1,K)*AS1(NX1,K)
     .                            -AS2(NX1,K)*PZ(NX1,K-1)*AS2(NX1,K)
     .                            -FX*AS1(2,K)*PZ(2,K)*AS1(2,K)      )
13    CONTINUE
      PZ(NX1,NY1) = 1./(AS3(NX1,NY1)
     .                      -AS1(NX1,NY1)*PZ(NX1-1,NY1)*AS1(NX1,NY1)
     .                      -AS2(NX1,NY1)*PZ(NX1,NY1-1)*AS2(NX1,NY1)
     .                      -FX*AS1(2,NY1)*PZ(2,NY1)*AS1(2,NY1)
     .                      -FY*AS2(NX1,2)*PZ(NX1,2)*AS2(NX1,2)      )
C
C----------------------------------------------------------------------
C.....boundary conditions:
C
      JA  =  1 + (1+LBOUN(1))/2 + (LBOUN(1)/3)*(NX-4)
      JE  = NX - (1+LBOUN(2))/2 - (LBOUN(2)/3)*(NX-4)
      KA  =  1 + (1+LBOUN(3))/2 + (LBOUN(3)/3)*(NY-4)
      KE  = NY - (1+LBOUN(4))/2 - (LBOUN(4)/3)*(NY-4)
      FXA = 1 - 2*(LBOUN(1)/2)+2*(LBOUN(1)/3)
      FXE = 1 - 2*(LBOUN(2)/2)+2*(LBOUN(2)/3)
      FYA = 1 - 2*(LBOUN(3)/2)+2*(LBOUN(3)/3)
      FYE = 1 - 2*(LBOUN(4)/2)+2*(LBOUN(4)/3)
C
C----------------------------------------------------------------------
C.....simple 5-point symmetric solution:
C
        IF(NS5 .EQ. 1) THEN
          CALL Z2SYM(ITMAX,ERROR,ITER,SER,LBOUN)
C
          DO 21 K=2,NY1
            F(1,K)  = FXA*F(JA,K)
            F(NX,K) = FXE*F(JE,K)
21        CONTINUE
          DO 22 J=2,NX1
            F(J,1)  = FYA*F(J,KA)
            F(J,NY) = FYE*F(J,KE)
22        CONTINUE
          F(1,1)    = FXA*FYA*F(JA,KA)
          F(1,NY)   = FXA*FYE*F(JA,KE)
          F(NX,1)   = FXE*FYA*F(JE,KA)
          F(NX,NY)  = FXE*FYE*F(JE,KE)
C
          RETURN
      ENDIF
C
C----------------------------------------------------------------------
C.....here solve general system with asymmetric terms and extra terms
C.....of 9-point system:
C
C.....store source term for general system:
      DO 31 K=2,NY1
      DO 31 J=2,NX1
        SA(J,K) = S(J,K)
31    CONTINUE
C
C.....convergence for general system:
      RS   = 0.0
      DO 32 K=2,NY1
      DO 32 J=2,NX1
        RS = RS + S(J,K)**2
32    CONTINUE
      ITER = 0
      ITT  = 1
      ERIN = ERROR
      SER0 = 0.1
C
C.....perform an iterative solution for skew-symmetric terms and extra
C.....terms of a 9-point stencil:
C
      DO 43 ITA=1,ITMAX
        CALL Z2ASY (JA,FXA,JE,FXE,KA,FYA,KE,FYE,ITT,LBOUN,
     .              A6,A7,A8,A9,NDX,NDY,AB1,s1,AB2,s2)
C.......Q = S-M*F
        CALL Z2MAT (Q,F,S,1,LBOUN)
        QS = 0.0
        DO 41 K=2,NY1
        DO 41 J=2,NX1
          QS = QS + Q(J,K)**2
41      CONTINUE
        SER = SQRT(QS/RS)
C
        IF(ITA .EQ. 1) THEN
          ERIN = AMAX1(ERROR,.5*AMIN1(SER,SER0))
        ELSE
          ERIN = AMAX1(ERROR,AMIN1(ERIN,SER**2/SER0))
        ENDIF
        SER0 = SER
C
        IF(LPRINT) WRITE (6,99) ITA,ITER,ITT,SER,ERIN
99      FORMAT (' Z2ASY: ITA,ITER,ITT,SER,ERIN',3I4,1P2E11.3)
C
        IF(SER .EQ. 0.0) RETURN
        IF ( ITER.NE.0.AND.(SER.LE.ERROR.OR.ITER.GE.ITMAX) ) RETURN
C
        IF ( LBOUN(1).GE.3.AND.LBOUN(3).GE.3 ) THEN
          DO 42 K=1,NY
          DO 42 J=1,NX
            F(J,K) = F(J,K) - F(2,2)
42        CONTINUE
        ENDIF
C
        CALL Z2SYM (ITMAX,ERIN,ITT,SER,LBOUN)
        ITER = ITER + ITT
43    CONTINUE
C
      RETURN
C
      END
C**********************************************************************
C**********************************************************************
C**********************************************************************
      SUBROUTINE Z2SYM (ITMAX,ERROR,ITER,SER,LBOUN)
      SAVE
C
C......dimension block
C
      REAL A6(NDX,NDY),A7(NDX,NDY),A8(NDX,NDY),A9(NDX,NDY)
        REAL AB1(NDY),AB2(NDY)
C
#include "z2param.inc"
C
*      DOUBLE PRECISION F,AS1,AA1,AS2,AA2,AS3,S,SA,PZ,R,P,Q,H,Z,G,    !VAX
*     .                 RS,RR1,PQ,A,QS,RR2,                          !VAX
*     .                 QQ,QP,Q1P,PP,W                               !VAX
C
      COMMON/ICCG1/ AS1(MX,MY),AA1(MX,MY),AS2(MX,MY),AA2(MX,MY),
     .              AS3(MX,MY),PZ(MX,MY),H(MX,MY),Z(MX,MY),G(MX,MY),
     .              NX,NY,NX1,NY1,LPRINT
      COMMON/ICCG2/ R(MX,MY),P(MX,MY),Q(MX,MY),
     .              F(MX,MY),S(MX,MY),SA(MX,MY)
C
      INTEGER LBOUN(4)
      LOGICAL LPRINT
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in z2code.F::Z2SYM(...)"
#endif
C----------------------------------------------------------------------
C
C.....R=S-M*F
      CALL Z2MAT (R,F,S,1,LBOUN)
C
C.....P=(L*D*LT)**-1*R
      CALL Z2INV (P,R,LBOUN)
C
      RS  = 0.0
      RR1 = 0.0
      DO 2 K=2,NY1
      DO 2 J=2,NX1
        RS  = RS  + S(J,K)*S(J,K)
        RR1 = RR1 + R(J,K)*P(J,K)
2     CONTINUE
C
C----------------------------------------------------------------------
      DO 16 ITER=1,ITMAX
C.......Q = M*P
        CALL Z2MAT (Q,P,S,2,LBOUN)
        PQ = 0.0
        DO 11 K=2,NY1
        DO 11 J=2,NX1
          PQ = PQ + Q(J,K)*P(J,K)
11      CONTINUE
        A = RR1/PQ
        DO 12 K=2,NY1
        DO 12 J=2,NX1
          F(J,K) = F(J,K) + A*P(J,K)
          R(J,K) = R(J,K) - A*Q(J,K)
12      CONTINUE
        IF ( MOD(ITER,2).EQ.0 ) THEN
C.........Q=S-M*F
          CALL Z2MAT (Q,F,S,1,LBOUN)
          QS = 0.0
          DO 13 K=2,NY1
          DO 13 J=2,NX1
            QS = QS + Q(J,K)**2
13        CONTINUE
          SER = SQRT(QS/RS)
          IF(LPRINT) WRITE(6,99) ITER,SER,ERROR
 99       FORMAT(' Z2SYM:        ITER,SER,ERROR',8X,I4,1P2E11.3)
          IF ( SER.LE.ERROR.OR.ITER.GE.ITMAX ) RETURN
        ENDIF
C.......Q = (L*D*LT)**-1*R
        CALL Z2INV (Q,R,LBOUN)
        RR2 = 0.0
        DO 14 K=2,NY1
        DO 14 J=2,NX1
          RR2 = RR2 + R(J,K)*Q(J,K)
14      CONTINUE
        A = RR2/RR1
        RR1 = RR2
        DO 15 K=2,NY1
        DO 15 J=2,NX1
          P(J,K) = Q(J,K) + A*P(J,K)
15      CONTINUE
16    CONTINUE
C
      RETURN
C
C**********************************************************************
      ENTRY Z2ASY (JA,FXA,JE,FXE,KA,FYA,KE,FYE,ITT,LBOUN,
     .             A6,A7,A8,A9,NDX,NDY,AB1,s1,AB2,s2)
C
        FB1 = s1
        FB2 = s2
        DO 20 K=2,NY1
          FB1 = FB1 + AB1(K)*F(2,K)
          FB2 = FB2 + AB2(K)*F(NX-1,K)
   20   CONTINUE
      DO 21 K=2,NY1
        F(1,K)  = FB1
        F(NX,K) = FB2
21    CONTINUE
      DO 22 J=2,NX1
        F(J,1)  = FYA*F(J,KA)
        F(J,NY) = FYE*F(J,KE)
22    CONTINUE
      F(1,1)    = FXA*FYA*F(JA,KA)
      F(1,NY)   = FXA*FYE*F(JA,KE)
      F(NX,1)   = FXE*FYA*F(JE,KA)
      F(NX,NY)  = FXE*FYE*F(JE,KE)
      DO 23 K = 1, NY
      DO 23 J = 1, NX
        G(J,K) = F(J,K)
23    CONTINUE
C
C----------------------------------------------------------------------
      DO 38 I= 1, ITT
C
        DO 31 K=2,NY1
        DO 31 J=2,NX1
          S(J,K) = SA(J,K)-AA1(J,K)*F(J-1,K)+AA1(J+1,K)*F(J+1,K)
     .                    -AA2(J,K)*F(J,K-1)+AA2(J,K+1)*F(J,K+1)
     .                    -A6(J,K)*F(J-1,K-1)-A7(J,K)*F(J-1,K+1)
     .                    -A8(J,K)*F(J+1,K-1)-A9(J,K)*F(J+1,K+1)
31      CONTINUE
C.......Q = S-M*F
        CALL Z2MAT (Q,F,S,1,LBOUN)
C.......R = (L*D*LT)**-1*Q
        CALL Z2INV (R,Q,LBOUN)
        RB1 = s1
        RB2 = s2
        DO 30 K=2,NY1
          RB1 = RB1 + AB1(K)*R(2,K)
          RB2 = RB2 + AB2(K)*R(NX-1,K)
   30   CONTINUE
        DO 32 K=2,NY1
          R(1,K)  = FXA*R(JA,K)
          R(NX,K) = FXE*R(JE,K)
32      CONTINUE
        DO 33 J=2,NX1
          R(J,1)  = FYA*R(J,KA)
          R(J,NY) = FYE*R(J,KE)
33      CONTINUE
        R(1,1)    = FXA*FYA*R(JA,KA)
        R(1,NY)   = FXA*FYE*R(JA,KE)
        R(NX,1)   = FXE*FYA*R(JE,KA)
        R(NX,NY)  = FXE*FYE*R(JE,KE)
C.......P = M*R
        CALL Z2MAT (P,R,S,2,LBOUN)
        DO 35 K = 2, NY1
        DO 35 J = 2, NX1
          P(J,K) = P(J,K)+AA1(J,K)*R(J-1,K)-AA1(J+1,K)*R(J+1,K)
     .                   +AA2(J,K)*R(J,K-1)-AA2(J,K+1)*R(J,K+1)
     .                   +A6(J,K)*R(J-1,K-1)+A7(J,K)*R(J-1,K+1)
     .                   +A8(J,K)*R(J+1,K-1)+A9(J,K)*R(J+1,K+1)
35      CONTINUE
        QQ  = 0.0
        QP  = 0.0
        Q1P = 0.0
        PP  = 0.0
        DO 36 K = 2, NY1
        DO 36 J = 2, NX1
          QQ = QQ + Q(J,K)*Q(J,K)
          QP  = QP + Q(J,K)*P(J,K)
          Q1P = Q1P + Z(J,K)*P(J,K)
          PP  = PP + P(J,K)*P(J,K)
36      CONTINUE
C
        A = QQ/QP
        IF(I .EQ. 1) THEN
          W = 1.0
        ELSE
c         DET = (1.+A*Q1P/QQ1)**2
c     .          -(1.-QQ/QQ1+2.*A*Q1P/QQ1+A*A*PP/QQ1)
c          IF(DET .GE. 0.0) THEN
c            W = ( 1.+A*(Q1P/QQ1)+SQRT(DET) )
c     .          /(1.-QQ/QQ1+2.*A*Q1P/QQ1+A*A*PP/QQ1)
c          ELSE
            W = (QQ1+A*Q1P)/(QQ1-QQ+2.*A*Q1P+A*A*PP)
c          ENDIF
        ENDIF
C
        DO 37 K=1,NY
        DO 37 J=1,NX
          Z(J,K) = Q(J,K)
          R(J,K) = W*A*R(J,K)+(W-1.0)*(F(J,K)-G(J,K))
          G(J,K) = F(J,K)
          F(J,K) = F(J,K)+R(J,K)
37      CONTINUE
C
        QQ1 = QQ
C
38    CONTINUE
C
C---------------------------------------------------------------------
      DO 51 K=2,NY1
      DO 51 J=2,NX1
          S(J,K) = SA(J,K)-AA1(J,K)*F(J-1,K)+AA1(J+1,K)*F(J+1,K)
     .                    -AA2(J,K)*F(J,K-1)+AA2(J,K+1)*F(J,K+1)
     .                    -A6(J,K)*F(J-1,K-1)-A7(J,K)*F(J-1,K+1)
     .                    -A8(J,K)*F(J+1,K-1)-A9(J,K)*F(J+1,K+1)
51    CONTINUE
C
      RETURN
C
      END










