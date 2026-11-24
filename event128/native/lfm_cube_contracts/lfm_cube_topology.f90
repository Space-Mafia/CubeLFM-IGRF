! Six-patch topology following the reviewed GAMERA-Cube signed frames.
! Coordinates and vector components are Cartesian; only indices are permuted.
module lfm_cube_topology
  use iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: cube_frame, cube_neighbor, cube_direction, cube_inverse, cube_vertex_key
  public :: cube_layout, cube_rank, cube_owner, cube_key_less, cube_entity_key
  integer, parameter :: frames(3,3,6)=reshape([ &
    1,0,0, 0,1,0, 0,0,1, &
   -1,0,0, 0,-1,0, 0,0,1, &
    0,1,0, -1,0,0, 0,0,1, &
    0,-1,0, 1,0,0, 0,0,1, &
    0,0,1, 1,0,0, 0,1,0, &
    0,0,-1, 1,0,0, 0,-1,0], [3,3,6])
  type cube_layout
    integer :: global(3)=0, dims(3)=0, coord(3)=0, first(3)=0, local(3)=0
    integer :: face=-1, rank=-1, size=0
  contains
    procedure :: init => layout_init
  end type
contains
  pure function cube_frame(face) result(q)
    integer, intent(in) :: face
    integer :: q(3,3)
    q=0
    if(face>=0.and.face<6) q=frames(:,:,face+1)
  end function

  pure subroutine cube_neighbor(face,edge,other,other_edge,reversed,status)
    integer, intent(in) :: face,edge ! 0:a-,1:a+,2:b-,3:b+
    integer, intent(out) :: other,other_edge,status
    logical, intent(out) :: reversed
    integer :: across,along,s,f,d,da,ds,orient,donor_along
    other=-1;other_edge=-1;reversed=.false.;status=1
    if(face<0.or.face>=6.or.edge<0.or.edge>=4) return
    across=2+edge/2;along=5-across;s=2*mod(edge,2)-1
    do f=0,5
      if(dot_product(frames(:,1,f+1),s*frames(:,across,face+1))==1) other=f
    enddo
    if(other<0.or.other==face) return
    da=0;ds=0
    do d=2,3
      orient=dot_product(frames(:,1,face+1),frames(:,d,other+1))
      if(orient/=0) then
        da=d;ds=orient
      endif
    enddo
    if(da==0) return
    other_edge=2*(da-2)+(ds+1)/2;donor_along=5-da
    reversed=dot_product(frames(:,along,face+1),frames(:,donor_along,other+1))<0
    status=0
  end subroutine

  pure subroutine cube_direction(face,a,b,x,status)
    integer,intent(in)::face
    real(real64),intent(in)::a,b
    real(real64),intent(out)::x(3)
    integer,intent(out)::status
    real(real64)::u,v,h
    x=0;status=1;h=acos(-1._real64)/4
    if(face<0.or.face>=6.or..not.all(ieee_is_finite([a,b]))) return
    if(max(abs(a),abs(b))>=2*h) return
    u=tan(a);v=tan(b)
    if(abs(a)==h) u=sign(1._real64,a)
    if(abs(b)==h) v=sign(1._real64,b)
    x=real(frames(:,1,face+1),real64)+u*frames(:,2,face+1)+v*frames(:,3,face+1)
    x=x/sqrt(sum(x*x));status=0
  end subroutine

  pure subroutine cube_inverse(x,face,a,b,status)
    real(real64),intent(in)::x(3)
    integer,intent(out)::face,status
    real(real64),intent(out)::a,b
    integer::f
    real(real64)::p,best,tol
    face=-1;status=1;a=0;b=0;best=0
    if(.not.all(ieee_is_finite(x))) return
    do f=0,5
      p=dot_product(x,real(frames(:,1,f+1),real64))
      tol=64*epsilon(p)*max(abs(p),abs(best))
      if(p>best+tol) then
        best=p;face=f
      endif
    enddo
    if(face<0) return
    a=atan2(dot_product(x,real(frames(:,2,face+1),real64)),best)
    b=atan2(dot_product(x,real(frames(:,3,face+1),real64)),best);status=0
  end subroutine

  pure subroutine cube_vertex_key(face,global,vertex,key,status)
    integer,intent(in)::face,global(3),vertex(3) ! zero based corner indices
    integer(int64),intent(out)::key(4)
    integer,intent(out)::status
    integer(int64)::n,a,b
    key=0;status=1
    if(face<0.or.face>=6.or.any(global<=0).or.global(2)/=global(3)) return
    if(any(vertex<0).or.any(vertex>global)) return
    n=int(global(2),int64);a=2_int64*vertex(2)-n;b=2_int64*vertex(3)-n
    key(1)=vertex(1)
    key(2:)=n*frames(:,1,face+1)+a*frames(:,2,face+1)+b*frames(:,3,face+1)
    status=0
  end subroutine

  pure logical function cube_key_less(a,b)
    integer(int64),intent(in)::a(:),b(:)
    integer::i
    cube_key_less=.false.
    do i=1,size(a)
      if(a(i)<b(i)) then
        cube_key_less=.true.;return
      else if(a(i)>b(i)) then
        return
      endif
    enddo
  end function

  pure subroutine cube_entity_key(vertices,key,orientation)
    ! Input is an oriented edge (2 vertices) or oriented face cycle (4).
    ! Canonical cycles preserve orientation, unlike an unordered vertex sort.
    integer(int64),intent(in)::vertices(:,:)
    integer(int64),intent(out)::key(17)
    integer,intent(out)::orientation
    integer::n,m,j,k,next,prev
    n=size(vertices,2);key=0;key(1)=n;orientation=1;m=1
    do j=2,n
      if(cube_key_less(vertices(:,j),vertices(:,m))) m=j
    enddo
    if(n==2) then
      if(m==2) orientation=-1
    else if(n==4) then
      next=mod(m,n)+1;prev=mod(m+n-2,n)+1
      if(cube_key_less(vertices(:,prev),vertices(:,next))) orientation=-1
    else
      orientation=0;return
    endif
    do j=0,n-1
      k=modulo(m-1+orientation*j,n)+1
      key(2+4*j:5+4*j)=vertices(:,k)
    enddo
  end subroutine

  subroutine layout_init(this,global,dims,rank,ranks,status)
    class(cube_layout),intent(out)::this
    integer,intent(in)::global(3),dims(3),rank,ranks
    integer,intent(out)::status
    integer::q,perface
    status=1
    if(any(dims<=0).or.any(global<=0)) return
    if(global(2)/=global(3).or.dims(2)/=dims(3)) return
    if(any(mod(global,dims)/=0)) return
    if(6_int64*product(int(dims,int64))/=int(ranks,int64)) return
    if(rank<0.or.rank>=ranks) return
    perface=product(dims);this%face=rank/perface;q=mod(rank,perface)
    this%coord(3)=mod(q,dims(3));q=q/dims(3)
    this%coord(2)=mod(q,dims(2));this%coord(1)=q/dims(2)
    this%global=global;this%dims=dims;this%local=global/dims
    this%first=this%coord*this%local;this%rank=rank;this%size=ranks;status=0
  end subroutine

  pure integer function cube_rank(face,coord,dims) result(rank)
    integer,intent(in)::face,coord(3),dims(3)
    rank=-1
    if(face<0.or.face>=6.or.any(dims<=0)) return
    if(any(coord<0).or.any(coord>=dims)) return
    rank=face*product(dims)+(coord(1)*dims(2)+coord(2))*dims(3)+coord(3)
  end function

  pure subroutine cube_owner(face,index,global,dims,rank,local,status)
    integer,intent(in)::face,index(3),global(3),dims(3)
    integer,intent(out)::rank,local(3),status
    integer::width(3)
    status=1;rank=-1;local=0
    if(any(dims<=0).or.any(global<=0)) return
    if(any(mod(global,dims)/=0).or.any(index<0).or.any(index>=global)) return
    width=global/dims
    rank=cube_rank(face,index/width,dims);local=mod(index,width)
    if(rank>=0) status=0
  end subroutine
end module
