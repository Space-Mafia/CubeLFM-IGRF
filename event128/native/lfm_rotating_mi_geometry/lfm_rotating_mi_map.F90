! Dynamic mapping adapters. SI, Cartesian vectors; no MHD/Poisson formula here.
module lfm_rotating_mi_map
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lfm_rotating_mi_geometry
  use lfm_tilt_geometry, only: dipole_field
  implicit none
  private
  integer,parameter::rk=real64
  real(rk),parameter::pi=acos(-1._rk)
  public::mi_cap_weights,mi_cap_evaluate,mi_dynamic_fac,mi_shell_scalar,mi_inner_potential
contains
  subroutine mi_cap_weights(x,y,p,ij,w,dw,status)
    ! Invert one native bilinear polar-plane element. dw is d shape/d(x,y).
    ! No nearest-point extrapolation. status2 means not covered by this mesh.
    real(rk),intent(in)::x(:,:),y(:,:),p(2)
    integer,intent(out)::ij(2),status
    real(rk),intent(out)::w(4),dw(2,4)
    real(rk)::v(2,4),a(2),b(2),c(2),z(2),ju(2),jv(2),det,u,t,du,dv,wu(4),wv(4)
    integer::j,k,it,l
    status=2;ij=0;w=0;dw=0
    if(any(shape(x)/=shape(y)).or..not.all(ieee_is_finite(p))) then
      status=1;return
    endif
    do k=1,size(x,2)-1
      do j=1,size(x,1)-1
        v(1,:)=[x(j,k),x(j+1,k),x(j,k+1),x(j+1,k+1)]
        v(2,:)=[y(j,k),y(j+1,k),y(j,k+1),y(j+1,k+1)]
        if(any(p<minval(v,dim=2)-1.e-13_rk).or.any(p>maxval(v,dim=2)+1.e-13_rk)) cycle
        a=v(:,2)-v(:,1);b=v(:,3)-v(:,1);c=v(:,4)-v(:,2)-v(:,3)+v(:,1)
        u=.5_rk;t=.5_rk
        do it=1,20
          z=v(:,1)+a*u+b*t+c*u*t-p
          ju=a+c*t;jv=b+c*u;det=ju(1)*jv(2)-ju(2)*jv(1)
          if(abs(det)<1.e-22_rk) exit
          du=(jv(2)*z(1)-jv(1)*z(2))/det
          dv=(-ju(2)*z(1)+ju(1)*z(2))/det
          u=u-du;t=t-dv
          if(max(abs(du),abs(dv))<1.e-12_rk) exit
          if(max(abs(u),abs(t))>100) exit
        enddo
        if(min(u,t)<-1.e-10_rk.or.max(u,t)>1+1.e-10_rk) cycle
        if(norm2(v(:,1)+a*u+b*t+c*u*t-p)>1.e-12_rk) cycle
        ju=a+c*t;jv=b+c*u;det=ju(1)*jv(2)-ju(2)*jv(1)
        if(abs(det)<1.e-22_rk) cycle
        w=[(1-u)*(1-t),u*(1-t),(1-u)*t,u*t]
        wu=[-(1-t),1-t,-t,t];wv=[-(1-u),-u,1-u,u]
        do l=1,4
          dw(:,l)=[jv(2)*wu(l)-ju(2)*wv(l),-jv(1)*wu(l)+ju(1)*wv(l)]/det
        enddo
        ij=[j,k];status=0;return
      enddo
    enddo
  end subroutine

  subroutine mi_cap_evaluate(r,q,rion,x,y,pn,ps,value,e,status)
    ! Phi is HELD in magnetic coordinates, not held on computational edges.
    ! Pull back the gradient of exactly the same bilinear Phi interpolation.
    real(rk),intent(in)::r(3),q(3,3),rion,x(:,:),y(:,:),pn(:,:),ps(:,:)
    real(rk),intent(out)::value,e(3)
    integer,intent(out)::status
    real(rk)::f(3),df(3,3),fm(3),rr,theta,cp,sp,w(4),dw(2,4),v(4),g(2),gm(3),et(3),ep(3)
    integer::h,ij(2),j,k
    value=0;e=0
    call mi_footpoint_jacobian(r,q,rion,f,df,h,status)
    if(status/=0) return
    fm=matmul(transpose(q),f);rr=norm2(fm(:2));theta=atan2(rr,abs(fm(3)))
    cp=1;sp=0
    if(rr>1.e-14_rk*rion) then
      cp=fm(1)/rr;sp=fm(2)/rr
    endif
    call mi_cap_weights(x,y,theta*[cp,sp],ij,w,dw,status)
    if(status/=0) return
    j=ij(1);k=ij(2)
    if(h>0) then
      v=[pn(j,k),pn(j+1,k),pn(j,k+1),pn(j+1,k+1)]
    else
      v=[ps(j,k),ps(j+1,k),ps(j,k+1),ps(j+1,k+1)]
    endif
    value=dot_product(w,v);g=matmul(dw,v)
    if(theta<1.e-12_rk) then
      gm=[g(1),g(2),0._rk]/rion
    else
      et=[cos(theta)*cp,cos(theta)*sp,-h*sin(theta)];ep=[-sp,cp,0._rk]
      gm=((g(1)*cp+g(2)*sp)*et+theta/sin(theta)*(-g(1)*sp+g(2)*cp)*ep)/rion
    endif
    e=mi_pullback_e(df,matmul(q,gm))
  end subroutine

  pure subroutine mi_shell_scalar(r,theta,phi,f,value,status,gradient,regular_ends)
    real(rk),intent(in)::r(3),theta(:),phi(:),f(:,:)
    real(rk),intent(out)::value
    integer,intent(out)::status
    real(rk),optional,intent(out)::gradient(3)
    logical,optional,intent(in)::regular_ends
    real(rk)::a,p,u,v,da,dp,et(3),ep(3),du,dv,s
    integer::j,k,n,m
    n=size(theta);m=size(phi);value=0;status=1
    if(present(gradient)) gradient=0
    if(n<2.or.m<3.or.any(shape(f)/=[n,m])) return
    if(any(theta(2:)<=theta(:n-1)).or.any(phi(2:)<=phi(:m-1))) return
    if(abs(phi(m)-phi(1)-2*pi)>1.e-10_rk) return
    a=atan2(norm2(r(2:)),r(1));p=modulo(atan2(r(3),r(2))-phi(1),2*pi)+phi(1)
    if(a<theta(1)-1.e-8_rk.or.a>theta(n)+1.e-8_rk) return
    j=1;k=1
    do while(j<n-1)
      if(theta(j+1)>=a) exit
      j=j+1
    enddo
    do while(k<m-1)
      if(phi(k+1)>=p) exit
      k=k+1
    enddo
    u=max(0._rk,min(1._rk,(a-theta(j))/(theta(j+1)-theta(j))))
    v=(p-phi(k))/(phi(k+1)-phi(k))
    du=1;dv=1
    if(present(regular_ends)) then
      if(regular_ends) then
        ! Minimal C1 endpoint Hermite reparameterization. All nodal Phi
        ! values are preserved. Only cells adjacent to the zero equator or
        ! zero polar-axis nodes change reconstruction. This enforces Phi=O(z^2)
        ! near magnetic equator, hence finite E_t/B_n without a B/velocity floor.
        if(j==1) then
          s=u;u=2*s*s-s*s*s;du=4*s-3*s*s
        else if(j==n-1) then
          s=1-u;u=1-(2*s*s-s*s*s);du=4*s-3*s*s
        endif
        if(abs(sin(phi(k)))<1.e-12_rk) then
          s=v;v=2*s*s-s*s*s;dv=4*s-3*s*s
        else if(abs(sin(phi(k+1)))<1.e-12_rk) then
          s=1-v;v=1-(2*s*s-s*s*s);dv=4*s-3*s*s
        endif
      endif
    endif
    value=(1-u)*(1-v)*f(j,k)+u*(1-v)*f(j+1,k)+(1-u)*v*f(j,k+1)+u*v*f(j+1,k+1)
    if(present(gradient)) then
      da=du*((1-v)*(f(j+1,k)-f(j,k))+v*(f(j+1,k+1)-f(j,k+1)))/(theta(j+1)-theta(j))
      dp=dv*((1-u)*(f(j,k+1)-f(j,k))+u*(f(j+1,k+1)-f(j+1,k)))/(phi(k+1)-phi(k))
      if(abs(sin(a))<1.e-12_rk) then
        if(da/=0.or.dp/=0) then
          status=2;return
        endif
      else
        et=[-sin(a),cos(a)*cos(p),cos(a)*sin(p)];ep=[0._rk,-sin(p),cos(p)]
        gradient=(da*et+dp/sin(a)*ep)/norm2(r)
      endif
    endif
    status=0
  end subroutine

  pure subroutine mi_inner_potential(r,q,rinner,theta,phi,table,value,e,status)
    ! Reconstruct the COMPLETE original INTRP inner-sphere table, including
    ! its zero-potential nodes, then rotate. Never hard-cut the cap interpolant.
    ! Along-dipole extension is used only for inner boundary and ghost points.
    ! status2: field line closes below the wall / exact equatorial limit;
    ! caller may use zero ONLY after validating the table's equatorial zero band.
    real(rk),intent(in)::r(3),q(3,3),rinner,theta(:),phi(:),table(:,:)
    real(rk),intent(out)::value,e(3)
    integer,intent(out)::status
    real(rk)::rm(3),sm(3),df(3,3),grad(3),dg(3),rr,r2,rho2,a,g,h,ratio
    integer::i,j
    value=0;e=0;status=1
    if(.not.all(ieee_is_finite(r)).or.rinner<=0) return
    rm=matmul(transpose(q),r);r2=sum(rm*rm);rr=sqrt(r2)
    if(rr<=0.or.rr>rinner*(1+1.e-10_rk)) return
    if(all(table==0)) then
      status=0;return
    endif
    rho2=sum(rm(:2)**2);ratio=rinner/rr
    ! On the exact spherical wall use its analytic identity, avoiding
    ! subtraction of two nearly equal O(1) terms as magnetic latitude ->0.
    if(abs(ratio-1)<1.e-13_rk) ratio=1
    g=(1-ratio)+ratio*(rm(3)/rr)**2
    if(g<=0.or.abs(rm(3))<32*epsilon(rr)*rr) then
      status=2;return
    endif
    h=sign(1._rk,rm(3));a=(rinner/rr)**1.5_rk
    sm=[a*rm(1),a*rm(2),h*rinner*sqrt(g)]
    call mi_shell_scalar(sm,theta,phi,table,value,status,grad,regular_ends=.true.)
    if(status/=0) return
    df=0
    do i=1,2
      do j=1,3
        df(i,j)=-1.5_rk*a*rm(i)*rm(j)/r2
      enddo
      df(i,i)=df(i,i)+a
    enddo
    dg=3*rinner*rho2*rm/rr**5;dg(:2)=dg(:2)-2*rinner*rm(:2)/rr**3
    df(3,:)=h*rinner*dg/(2*sqrt(g))
    e=-matmul(q,matmul(transpose(df),grad))
    status=0
  end subroutine

  subroutine mi_dynamic_fac(x,y,q,rion,mu,theta,radii,phi,jcart,smooth,rn,rs,bad,status)
    ! Jcart: native JPARA Cartesian current in A/m^2, including periodic node.
    ! First project J.B/B^2 on native shell; optional donor 9-point smoothing
    ! retains the original logical interior extent and unchanged boundary rows.
    ! Then sample this invariant at each magnetic-cap node's physical source.
    ! RHS is signed OUTWARD radial current at ionospheric foot (both caps).
    real(rk),intent(in)::x(:,:),y(:,:),q(3,3),rion,mu,theta(:),radii(:),phi(:),jcart(:,:,:)
    logical,intent(in)::smooth
    real(rk),intent(out)::rn(:,:),rs(:,:)
    integer,intent(out)::bad,status
    real(rk)::a(size(theta),size(phi)),ac(size(theta),size(phi)),r(3),b(3),f(3),fm(3),t,p,z
    integer::j,k,h,st,n,m
    rn=0;rs=0;bad=0;status=1;n=size(theta);m=size(phi)
    if(any(shape(jcart)/=[3,n,m]).or.size(radii)/=n) return
    if(any(shape(x)/=shape(y)).or.any(shape(rn)/=shape(x)).or.any(shape(rs)/=shape(x))) return
    if(.not.all(ieee_is_finite(jcart)).or.mu==0) return
    do k=1,m
      do j=1,n
        r=radii(j)*[cos(theta(j)),sin(theta(j))*cos(phi(k)),sin(theta(j))*sin(phi(k))]
        b=dipole_field(r,mu*q(:,3))
        a(j,k)=dot_product(jcart(:,j,k),b)/sum(b*b)
      enddo
    enddo
    if(smooth) then
      ac=a
      do k=2,m-2
        do j=2,n-1
          ac(j,k)=.25_rk*a(j,k)+.125_rk*(a(j+1,k)+a(j,k+1)+a(j-1,k)+a(j,k-1)) &
            +.0625_rk*(a(j-1,k-1)+a(j-1,k+1)+a(j+1,k-1)+a(j+1,k+1))
        enddo
      enddo
      a=ac
    endif
    do h=-1,1,2
      do k=3,size(x,2)-2
        do j=2,size(x,1)-1
          t=hypot(x(j,k),y(j,k));p=atan2(y(j,k),x(j,k))
          fm=rion*[sin(t)*cos(p),sin(t)*sin(p),h*cos(t)];f=matmul(q,fm)
          call mi_shell_intersection(f,q,rion,theta,radii,r,st)
          if(st/=0) then
            bad=bad+1;cycle
          endif
          call mi_shell_scalar(r,theta,phi,a,z,st)
          if(st/=0) then
            bad=bad+1;cycle
          endif
          b=dipole_field(f,mu*q(:,3));z=z*dot_product(b,f/rion)
          if(h>0) rn(j,k)=z
          if(h<0) rs(j,k)=z
        enddo
      enddo
    enddo
    status=0
    if(bad>0) status=3 ! Never silently lose an unmapped cap current.
  end subroutine
end module
