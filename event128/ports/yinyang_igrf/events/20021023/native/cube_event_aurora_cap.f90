module cube_event_aurora_cap
 use iso_c_binding
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 use lfm_cube_yy_native_cap_frame,only:yy_native_cap_frame
 use lfm_cube_exact_cap_lookup,only:exact_cap_configure,exact_cap_weights
 use lfm_rotating_mi_map,only:mi_cap_weights
 implicit none
 private
 public::event_cap_update, event_np, event_channels
 ! Sample DPB evidence on the released two-cap reference MLAT/MLT grid.
 integer,parameter::nt=99,event_np=384,event_channels=12
 type,bind(C)::auroral_point
  real(c_double)::ep,eh,ap,ah,number,energy,mean,sp,sh
  integer(c_int)::diffuse
 end type
 interface
  subroutine feet(nt,np,t,v) bind(C,name='cube_event_dpb_feet')
   import
   integer(c_int),value::nt,np
   real(c_double),value::t
   real(c_double),intent(out)::v(*)
  end subroutine
  integer(c_int) function advance_dpb(nt,np,f,p,ready,r,edge,mask) bind(C,name='cube_event_dpb_step')
   import
   integer(c_int),value::nt,np
   integer(c_int),intent(inout)::ready(*)
   real(c_double),intent(in)::f(*),p(*)
   real(c_double),intent(inout)::r(*),edge(*)
   real(c_double),intent(out)::mask(*)
  end function
  integer(c_int) function masks(n,np,t,v,edge,m,z,lat,mlt) bind(C,name='cube_event_geographic_mask')
   import
   integer(c_int),value::n,np
   real(c_double),value::t
   real(c_double),intent(in)::v(*),edge(*)
   real(c_double),intent(out)::m(*),z(*),lat(*),mlt(*)
  end function
  integer(c_int) function precip(t,z,valid,rho,cs,fac,inc,m,ready,oldp,oldh,pt) bind(C,name='cube_yy_aurora_point')
   import
   integer(c_int),value::valid,ready
   real(c_double),value::t,z,rho,cs,fac,inc,m,oldp,oldh
   type(auroral_point),intent(out)::pt
  end function
 end interface
contains
 subroutine event_cap_update(t,x,y,f,p,sources,trace,ready,radius,edge,sp,sh,diag,status)
  real(c_double),intent(in)::t,x(:,:),y(:,:),f(:,:,:),p(:,:,:),sources(:,:,:,:)
  integer,intent(in)::trace(:,:,:)
  integer(c_int),intent(inout)::ready(2)
  real(c_double),intent(inout)::radius(2),edge(event_np,2),sp(:,:,:),sh(:,:,:)
  real(c_double),intent(out)::diag(:,:,:,:)
  integer,intent(out)::status
  real(c_double),allocatable::points(:,:,:,:),fac(:,:,:),pot(:,:,:),mask(:,:,:),direction(:,:,:,:),cut(:,:,:),zen(:,:,:)
  real(c_double),allocatable::lat(:,:,:),mlt(:,:,:)
  real(c_double)::frame(3,3),v(3),u(2),w(4),dw(2,4),a,b,theta,phi,pi
  integer::i,j,h,k,l,side,ij(2),st,was_ready(2)
  logical::exact
  type(auroral_point)::pt
  status=1; pi=acos(-1.d0);was_ready=ready
  frame=yy_native_cap_frame(t)
  allocate(points(3,event_np,nt,2),fac(event_np,nt,2),pot(event_np,nt,2),mask(event_np,nt,2))
  call feet(nt,event_np,t,points)
  call exact_cap_configure(x,y,exact)
  do h=1,2
   do i=1,nt
    do j=1,event_np
     v=matmul(transpose(frame),points(:,j,i,h))
     theta=atan2(norm2(v(:2)),abs(v(3)));phi=atan2(v(2),v(1))
     u=theta*[cos(phi),sin(phi)]
     if(exact)then
      call exact_cap_weights(u,ij,w,dw,st)
     else
      call mi_cap_weights(x,y,u,ij,w,dw,st)
     endif
     if(st==2)then
      fac(j,i,h)=0;pot(j,i,h)=0
     else if(st==0)then
      k=ij(1);l=ij(2);side=merge(1,2,v(3)>0)
      fac(j,i,h)=dot_product(w,[f(k,l,side),f(k+1,l,side),f(k,l+1,side),f(k+1,l+1,side)])
      pot(j,i,h)=dot_product(w,[p(k,l,side),p(k+1,l,side),p(k,l+1,side),p(k+1,l+1,side)])
     else
      return
     endif
    enddo
   enddo
  enddo
  if(advance_dpb(nt,event_np,fac,pot,ready,radius,edge,mask)/=0)return
  allocate(direction(3,size(x,1),size(x,2),2),cut(size(x,1),size(x,2),2),zen(size(x,1),size(x,2),2))
  do h=1,2
   do j=1,size(x,2)
    do i=1,size(x,1)
     theta=hypot(x(i,j),y(i,j));phi=atan2(y(i,j),x(i,j))
     if(h==2) theta=pi-theta
     direction(:,i,j,h)=matmul(frame,[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)])
    enddo
   enddo
  enddo
  allocate(lat(size(x,1),size(x,2),2),mlt(size(x,1),size(x,2),2))
  if(masks(2*size(x),event_np,t,direction,edge,cut,zen,lat,mlt)/=0)return
  do h=1,2
   do j=1,size(x,2)
    do i=1,size(x,1)
     a=sp(i,j,h);b=sh(i,j,h)
     st=precip(t,zen(i,j,h),merge(1,0,trace(i,j,h)==0),sources(1,i,j,h),sources(2,i,j,h), &
        sources(3,i,j,h),sources(4,i,j,h),cut(i,j,h),was_ready(h),a,b,pt)
     if(st/=0)return
     sp(i,j,h)=pt%sp;sh(i,j,h)=pt%sh
     diag(:,i,j,h)=[pt%energy,pt%number,pt%mean,pt%ep,pt%eh,pt%ap,pt%ah,cut(i,j,h),zen(i,j,h),real(pt%diffuse,c_double), &
                   lat(i,j,h),mlt(i,j,h)]
    enddo
   enddo
  enddo
  if(.not.all(ieee_is_finite(diag)))return
  status=0
 end subroutine
end module
