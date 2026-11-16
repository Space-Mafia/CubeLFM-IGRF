      subroutine cube_yy_g2int(fnorm,FPRESS,fx,fy,fz,fnx,fny,fnz,face,
     $     x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3)
      SAVE
C     This subroutine does Gaussian integration of the function 
C
C     The Gaussian integration uses the first twelve Legendre polynomials
C     as the basis functions. Abromowitz and Stegun page 916.
#ifndef T3E
      double precision bx,by,bz,bsq
#endif
      real dx(3),dy(3),dz(3)
      real WT(12),A(12)
      DATA A/0.1252334085,  0.3678314989,  0.5873179542,
     $     0.7699026741,  0.9041172563,  0.9815606342,
     $     -0.1252334085, -0.3678314989, -0.5873179542,
     $     -0.7699026741, -0.9041172563, -0.9815606342/
      DATA WT/0.2491470458,  0.2334925365,  0.2031674267,
     $     0.1600783285,  0.1069393259,  0.0471753363,
     $     0.2491470458, 0.2334925365, 0.2031674267,
     $     0.1600783285, 0.1069393259,  0.0471753363/
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subdipol.F::cube_yy_g2int(...)"
#endif
      dx(1) = x1-x0
      dx(2) = x2-x0
      dx(3) = x3+x0-x2-x1
      dy(1) = y1-y0
      dy(2) = y2-y0
      dy(3) = y3+y0-y2-y1
      dz(1) = z1-z0
      dz(2) = z2-z0
      dz(3) = z3+z0-z2-z1
C
      fnorm = 0.
      FPRESS = 0.
      fx = 0.
      fy = 0.
      fz = 0.
      fnx = 0.
      fny = 0.
      fnz = 0.
C     
      do 200 i=1,12
         do 200 j=1,12
            eta = 0.5*(1. + a(i))
            psi = 0.5*(1. + a(j))
C
            x = x0 + dx(1)*eta + dx(2)*psi + dx(3)*eta*psi
            y = y0 + dy(1)*eta + dy(2)*psi + dy(3)*eta*psi
            z = z0 + dz(1)*eta + dz(2)*psi + dz(3)*eta*psi
C
            xeta = dx(1) + dx(3)*psi
            yeta = dy(1) + dy(3)*psi
            zeta = dz(1) + dz(3)*psi
C     
            xpsi = dx(2) + dx(3)*eta
            ypsi = dy(2) + dy(3)*eta
            zpsi = dz(2) + dz(3)*eta
C
C  note the fancy footwork with 1.e10 to make the vaxvms version behave
C
            xn = 1.e-10*(yeta*zpsi-zeta*ypsi)
            yn = 1.e-10*(zeta*xpsi-xeta*zpsi)
            zn = 1.e-10*(xeta*ypsi-xpsi*yeta)
            bx = cube_yy_bxqq0(x,y,z)
            by = cube_yy_byqq0(x,y,z)
            bz = cube_yy_bzqq0(x,y,z)
C
            bdotn = 2.5e9*(xn*bx + yn*by + zn*bz)*wt(i)*wt(j)
            rnsq = xn*xn+yn*yn+zn*zn
            rn = sqrt(rnsq)
            area = 2.5e9*rn*wt(i)*wt(j)
            BSQ = (BX*BX+BY*BY+BZ*BZ)
            fnorm = fnorm + bdotn
            fpress = fpress + area*BSQ
            fx = fx + area*bx
            fy = fy + area*by
            fz = fz + area*bz
            fnx = fnx + bdotn*bx
            fny = fny + bdotn*by
            fnz = fnz + bdotn*bz
C
  200   continue
C
        fnx = fnx/(fnorm + 1.e-20)
        fny = fny/(fnorm + 1.e-20)
        fnz = fnz/(fnorm + 1.e-20)
        fnorm = fnorm/face
        fpress = fpress/face
        fx = fx/face
        fy = fy/face
        fz = fz/face
C
C
        return
        end
C
