! Exact installed-geometry index for the existing native bilinear cell.
! No new interpolation, extrapolation, mapping, solver, cap or boundary.
module lfm_cube_exact_cap_lookup
 use iso_fortran_env,only:real64,int64
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_round_type,ieee_get_rounding_mode,operator(==)
 use lfm_rotating_mi_map,only:original_cap_weights=>mi_cap_weights
 implicit none
 private
 integer,parameter::rk=real64,index_side=64,index_count=index_side*index_side
 type cap_index
  logical::ready=.false.
  type(ieee_round_type)::rounding
  real(rk),allocatable::x(:,:),y(:,:),v(:,:,:),low(:,:),high(:,:)
  integer,allocatable::j(:),k(:),offset(:),cell(:)
  real(rk)::origin(2),scale(2),upper(2)
 end type
 type(cap_index),save::active
 public::exact_cap_configure,exact_cap_weights
contains
 pure integer function bucket(value,origin,scale)result(b)
  real(rk),intent(in)::value,origin,scale
  b=max(0,min(index_side-1,floor((value-origin)*scale)))+1
 end function
 subroutine exact_cap_configure(x,y,usable)
  real(rk),intent(in)::x(:,:),y(:,:)
  logical,intent(out)::usable
  type(cap_index)::candidate
  type(ieee_round_type)::mode
  integer::nc,nj,nk,j,k,c,i,bx,by,lowbin(2),highbin(2),b,st
  integer::counts(index_count),cursor(index_count)
  integer(int64)::total
  real(rk)::span(2)
  usable=.false.
  if(any(shape(x)/=shape(y)))return
  if(min(size(x,1),size(x,2))<2)return
  if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(y)))return
  call ieee_get_rounding_mode(mode)
  if(active%ready)then
   if(all(shape(x)==shape(active%x)))then
    if(mode==active%rounding)then
     if(all(transfer(x,[0_int64],size(x))==transfer(active%x,[0_int64],size(x))).and. &
        all(transfer(y,[0_int64],size(y))==transfer(active%y,[0_int64],size(y))))then
      usable=.true.;return
     endif
    endif
   endif
  endif
  active%ready=.false.
  ! Cache eligibility only: extreme coordinates use the original routine.
  ! No input is rejected and no new physical boundary is introduced.
  if(maxval(abs(x))>huge(1._rk)/16.or.maxval(abs(y))>huge(1._rk)/16)return
  nj=size(x,1);nk=size(x,2);nc=(nj-1)*(nk-1)
  allocate(candidate%x(nj,nk),candidate%y(nj,nk),candidate%v(2,4,nc), &
   candidate%low(2,nc),candidate%high(2,nc),candidate%j(nc),candidate%k(nc), &
   candidate%offset(index_count+1),stat=st)
  if(st/=0)return
  candidate%x=x;candidate%y=y;candidate%rounding=mode;c=0
  do k=1,nk-1
   do j=1,nj-1
    c=c+1;candidate%j(c)=j;candidate%k(c)=k
    candidate%v(1,:,c)=[x(j,k),x(j+1,k),x(j,k+1),x(j+1,k+1)]
    candidate%v(2,:,c)=[y(j,k),y(j+1,k),y(j,k+1),y(j+1,k+1)]
    candidate%low(:,c)=minval(candidate%v(:,:,c),dim=2)-1.e-13_rk
    candidate%high(:,c)=maxval(candidate%v(:,:,c),dim=2)+1.e-13_rk
   enddo
  enddo
  candidate%origin=minval(candidate%low,dim=2);candidate%upper=maxval(candidate%high,dim=2)
  span=candidate%upper-candidate%origin
  if(any(span<=0).or..not.all(ieee_is_finite(span)))return
  candidate%scale=index_side/span
  counts=0
  do c=1,nc
   do i=1,2
    lowbin(i)=bucket(candidate%low(i,c),candidate%origin(i),candidate%scale(i))
    highbin(i)=bucket(candidate%high(i,c),candidate%origin(i),candidate%scale(i))
   enddo
   do by=lowbin(2),highbin(2);do bx=lowbin(1),highbin(1)
    b=bx+index_side*(by-1);counts(b)=counts(b)+1
   enddo;enddo
  enddo
  total=sum(int(counts,int64));if(total>huge(1)-1)return
  allocate(candidate%cell(int(total)),stat=st);if(st/=0)return
  candidate%offset(1)=1
  do b=1,index_count
   candidate%offset(b+1)=candidate%offset(b)+counts(b)
  enddo
  cursor=candidate%offset(:index_count)
  ! Insert in original k/j order, so shared-edge ties choose the same cell.
  do c=1,nc
   do i=1,2
    lowbin(i)=bucket(candidate%low(i,c),candidate%origin(i),candidate%scale(i))
    highbin(i)=bucket(candidate%high(i,c),candidate%origin(i),candidate%scale(i))
   enddo
   do by=lowbin(2),highbin(2);do bx=lowbin(1),highbin(1)
    b=bx+index_side*(by-1);candidate%cell(cursor(b))=c;cursor(b)=cursor(b)+1
   enddo;enddo
  enddo
  call move_alloc(candidate%x,active%x);call move_alloc(candidate%y,active%y)
  call move_alloc(candidate%v,active%v);call move_alloc(candidate%low,active%low);call move_alloc(candidate%high,active%high)
  call move_alloc(candidate%j,active%j);call move_alloc(candidate%k,active%k)
  call move_alloc(candidate%offset,active%offset);call move_alloc(candidate%cell,active%cell)
  active%rounding=mode;active%origin=candidate%origin;active%upper=candidate%upper;active%scale=candidate%scale
  active%ready=.true.;usable=.true.
 end subroutine
 subroutine exact_cap_weights(p,ij,w,dw,status)
  real(rk),intent(in)::p(2)
  integer,intent(out)::ij(2),status
  real(rk),intent(out)::w(4),dw(2,4)
  real(rk)::v(2,4),a(2),b(2),c(2),z(2),ju(2),jv(2),det,u,t,du,dv,wu(4),wv(4)
  integer::j,k,it,l,bin,entry,cell
  type(ieee_round_type)::mode
  status=2;ij=0;w=0;dw=0
  if(.not.all(ieee_is_finite(p)))then
   status=1;return
  endif
  if(.not.active%ready)then
   status=1;return
  endif
  call ieee_get_rounding_mode(mode)
  if(.not.(mode==active%rounding))then
   call original_cap_weights(active%x,active%y,p,ij,w,dw,status);return
  endif
  if(any(p<active%origin).or.any(p>active%upper))return
  bin=bucket(p(1),active%origin(1),active%scale(1))+index_side*(bucket(p(2),active%origin(2),active%scale(2))-1)
  do entry=active%offset(bin),active%offset(bin+1)-1
   cell=active%cell(entry)
   if(any(p<active%low(:,cell)).or.any(p>active%high(:,cell)))cycle
   v=active%v(:,:,cell);j=active%j(cell);k=active%k(cell)
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
 end subroutine
end module
