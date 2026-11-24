! Native hexahedral geometry and integral CT contracts, independent of MHD.
module lfm_cube_mesh
  use iso_fortran_env, only: real64,int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lfm_cube_topology
  use lfm_cube_shell_map
  implicit none
  private
  public :: cube_patch, cube_meridian, cube_cross, cube_quad, cube_curl_update, cube_divergence, cube_bzz_calibration
  public :: cube_cell_statistics, cube_inverse3, cube_polynomial_moments
  type cube_meridian
    ! Small axisymmetric source table only. Never replicate a 3D state per rank.
    real(real64),allocatable :: xy(:,:,:),theta(:,:)
    integer :: ni=0,nj=0
    type(lfm_shell_map)::shell
    real(real64)::physical_scale=1
  contains
    procedure :: init => meridian_init
    procedure :: load => meridian_load
    procedure :: point => meridian_point
    procedure :: cell_vertices => meridian_cell_vertices
    procedure :: scale => meridian_scale
    procedure :: reference => meridian_reference
  end type
  type cube_patch
    type(cube_layout) :: layout
    real(real64),allocatable :: x(:,:,:,:),area(:,:,:,:,:),volume(:,:,:)
    real(real64),allocatable :: legacy_volume(:,:,:),legacy_length(:,:,:),cfl_length(:,:,:)
  contains
    procedure :: init => patch_init
    procedure :: init_mapped => patch_init_mapped
    procedure :: face_vertices
    procedure :: entity => patch_entity
    procedure :: edge_ends
  end type
contains
  pure function cube_cross(a,b) result(c)
    real(real64),intent(in)::a(3),b(3)
    real(real64)::c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
  pure function cube_quad(v) result(a)
    real(real64),intent(in)::v(3,4)
    real(real64)::a(3)
    a=.5_real64*cube_cross(v(:,3)-v(:,1),v(:,4)-v(:,2))
  end function
  pure real(real64) function volume_moment(v) result(w)
    real(real64),intent(in)::v(3,4)
    real(real64)::c(3)
    integer::j,k
    c=sum(v,dim=2)/4;w=0
    do j=1,4
      k=mod(j,4)+1
      w=w+dot_product(c,cube_cross(v(:,j),v(:,k)))/6
    enddo
  end function
  subroutine meridian_init(m,xy,status)
    class(cube_meridian),intent(out)::m
    real(real64),intent(in)::xy(:,0:,0:)
    integer,intent(out)::status
    integer::i
    status=1
    if(size(xy,1)/=2.or.size(xy,2)<2.or.size(xy,3)<3) return
    if(.not.all(ieee_is_finite(xy)).or.any(xy(2,:,:)<0)) return
    m%ni=size(xy,2)-1;m%nj=size(xy,3)-1
    if(any(xy(2,:,0)/=0).or.any(xy(2,:,m%nj)/=0)) return
    if(any(xy(1,:,0)<=0).or.any(xy(1,:,m%nj)>=0)) return
    allocate(m%xy(2,0:m%ni,0:m%nj),m%theta(0:m%ni,0:m%nj))
    m%xy=xy;m%theta=atan2(xy(2,:,:),xy(1,:,:))
    do i=0,m%ni
      if(any(m%theta(i,1:)<=m%theta(i,:m%nj-1))) return
    enddo
    status=0
  end subroutine
  subroutine meridian_load(m,path,status,radial_cells)
    class(cube_meridian),intent(out)::m
    character(*),intent(in)::path
    integer,intent(out)::status
    integer,intent(in),optional::radial_cells
    real(real64),allocatable::xy(:,:,:)
    integer::u,st,ni,nj,i,j
    character(80)::magic,extra
    real(real64)::direction(3),point(3),angle
    status=1
    open(newunit=u,file=path,status='old',action='read',iostat=st)
    if(st/=0) return
    read(u,'(a)',iostat=st) magic
    if(st==0.and.(trim(magic)=='GAMERA_CS_LFMD_V1'.or.trim(magic)=='LFM_CUBE_SHELL_V1')) then
      close(u);call m%shell%load(path,st)
      if(st/=0) return
      m%ni=m%shell%ns-1;m%nj=(m%shell%nt-1)/2
      if(present(radial_cells)) m%ni=radial_cells
      if(m%ni<1.or.mod(m%shell%nt-1,2)/=0) return
      allocate(m%xy(2,0:m%ni,0:m%nj),m%theta(0:m%ni,0:m%nj))
      do i=0,m%ni
        do j=0,m%nj
          angle=m%shell%theta(j);direction=[cos(angle),sin(angle),0._real64]
          call m%shell%forward(real(i,real64)/m%ni,direction,point,st)
          if(st/=0) return
          m%xy(:,i,j)=point(:2);m%theta(i,j)=atan2(point(2),point(1))
        enddo
      enddo
      status=0;return
    endif
    if(st==0.and.trim(magic)=='LFM_CUBE_MERIDIAN_V1') then
      read(u,*,iostat=st) ni,nj
      if(st==0) then
        if(ni>0.and.ni<=100000.and.nj>1.and.nj<=10000) then
          allocate(xy(2,0:ni,0:nj))
          do i=0,ni
            do j=0,nj
              read(u,*,iostat=st) xy(:,i,j)
              if(st/=0) exit
            enddo
            if(st/=0) exit
          enddo
          if(st==0) then
            read(u,'(a)',iostat=st) extra
            if(st<0) call m%init(xy,status)
          endif
        endif
      endif
    endif
    close(u)
  end subroutine
  subroutine meridian_point(m,face,nang,vertex,x,status)
    class(cube_meridian),intent(in)::m
    integer,intent(in)::face,nang,vertex(3)
    real(real64),intent(out)::x(3)
    integer,intent(out)::status
    real(real64)::d(3),a,b,t,px(2),edge(2),den,r,pi
    integer::i,lo,hi,mid,st
    status=1;x=0;i=vertex(1);pi=acos(-1._real64)
    if(.not.allocated(m%xy).or.nang<=0) return
    a=(2._real64*vertex(2)/nang-1)*pi/4
    b=(2._real64*vertex(3)/nang-1)*pi/4
    call cube_direction(face,a,b,d,st)
    if(st/=0) return
    if(m%shell%ns>0) then
      call m%shell%forward(real(i,real64)/m%ni,d,x,status);x=x*m%physical_scale;return
    endif
    if(i<0.or.i>m%ni) return
    t=atan2(sqrt(sum(d(2:)**2)),d(1));lo=0;hi=m%nj
    do while(hi-lo>1)
      mid=(lo+hi)/2
      if(t>=m%theta(i,mid)) then
        lo=mid
      else
        hi=mid
      endif
    enddo
    px=m%xy(:,i,lo);edge=m%xy(:,i,lo+1)-px
    den=cos(t)*edge(2)-sin(t)*edge(1)
    if(abs(den)<=tiny(den)) return
    r=(px(1)*edge(2)-px(2)*edge(1))/den
    if(.not.ieee_is_finite(r).or.r<=0) return
    x=r*d;status=0
  end subroutine
  subroutine meridian_scale(m,factor)
    class(cube_meridian),intent(inout)::m
    real(real64),intent(in)::factor
    if(.not.ieee_is_finite(factor).or.factor<=0) error stop 'positive geometry scale required'
    m%xy=m%xy*factor;m%physical_scale=m%physical_scale*factor
  end subroutine
  subroutine meridian_reference(m,x,t,direction,status)
    class(cube_meridian),intent(in)::m
    real(real64),intent(in)::x(3)
    real(real64),intent(out)::t,direction(3)
    integer,intent(out)::status
    if(m%shell%ns>0) then
      call m%shell%inverse(x/m%physical_scale,t,direction,status)
    else
      t=0;direction=0;status=1
      if(norm2(x)<=0.or..not.all(ieee_is_finite(x))) return
      direction=x/norm2(x);status=0
    endif
  end subroutine
  subroutine meridian_cell_vertices(m,face,nang,index,v,status)
    class(cube_meridian),intent(in)::m
    integer,intent(in)::face,nang,index(3) ! one-based cells; angular ghosts allowed
    real(real64),intent(out)::v(3,0:1,0:1,0:1)
    integer,intent(out)::status
    integer::a,b,c,st
    status=1
    do c=0,1
      do b=0,1
        do a=0,1
          call m%point(face,nang,index-1+[a,b,c],v(:,a,b,c),st)
          if(st/=0) return
        enddo
      enddo
    enddo
    status=0
  end subroutine
  subroutine patch_init_mapped(p,layout,meridian,status)
    class(cube_patch),intent(out)::p
    type(cube_layout),intent(in)::layout
    type(cube_meridian),intent(in)::meridian
    integer,intent(out)::status
    ! The log-shell radii are unused in the explicitly selected mapped branch.
    call patch_init(p,layout,1._real64,2._real64,status,meridian)
  end subroutine
  subroutine patch_init(p,layout,rin,rout,status,meridian)
    class(cube_patch),intent(out)::p
    type(cube_layout),intent(in)::layout
    real(real64),intent(in)::rin,rout
    integer,intent(out)::status
    type(cube_meridian),optional,intent(in)::meridian
    integer::i,j,k,d,n(3),g(3),st
    real(real64)::a,b,r,v(3,4),pos(3),pi,delta(3,3),centre(3),face_max
    status=1
    if(.not.all(ieee_is_finite([rin,rout]))) return
    if(rin<=0.or.rout<=rin.or.layout%face<0) return
    if(present(meridian)) then
      if(.not.allocated(meridian%xy).or.layout%global(1)/=meridian%ni) return
    endif
    p%layout=layout;n=layout%local;pi=acos(-1._real64)
    allocate(p%x(3,0:n(1),0:n(2),0:n(3)),p%area(3,0:n(1),0:n(2),0:n(3),3))
    allocate(p%volume(n(1),n(2),n(3)),p%legacy_volume(n(1),n(2),n(3)), &
      p%legacy_length(n(1),n(2),n(3)),p%cfl_length(n(1),n(2),n(3)));p%area=0
    do k=0,n(3)
      do j=0,n(2)
        do i=0,n(1)
          g=layout%first+[i,j,k]
          if(present(meridian)) then
            call meridian%point(layout%face,layout%global(2),g,pos,st)
            if(st/=0) return
            p%x(:,i,j,k)=pos
          else
            a=(2._real64*g(2)/layout%global(2)-1)*pi/4
            b=(2._real64*g(3)/layout%global(3)-1)*pi/4
            call cube_direction(layout%face,a,b,pos,st)
            if(st/=0) return
            r=rin*exp(log(rout/rin)*real(g(1),real64)/layout%global(1))
            p%x(:,i,j,k)=r*pos
          endif
        enddo
      enddo
    enddo
    do d=1,3
      do k=0,n(3)
        do j=0,n(2)
          do i=0,n(1)
            if(d/=1.and.i==0) cycle
            if(d/=2.and.j==0) cycle
            if(d/=3.and.k==0) cycle
            call p%face_vertices(d,[i,j,k],v)
            p%area(:,i,j,k,d)=cube_quad(v)
          enddo
        enddo
      enddo
    enddo
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p%volume(i,j,k)=0
          do d=1,3
            g=[i,j,k];call p%face_vertices(d,g,v)
            centre=sum(v,dim=2)/4
            p%volume(i,j,k)=p%volume(i,j,k)+volume_moment(v)
            g(d)=g(d)-1;call p%face_vertices(d,g,v)
            delta(:,d)=centre-sum(v,dim=2)/4
            p%volume(i,j,k)=p%volume(i,j,k)-volume_moment(v)
          enddo
          p%legacy_volume(i,j,k)=abs(dot_product(delta(:,1),cube_cross(delta(:,2),delta(:,3))))
          ! Preserve metric.F's five-face CFL expression, including lower-K
          ! only. A symmetric six-face maximum would be a separate algorithm.
          face_max=max(norm2(p%area(:,i,j,k,1)),norm2(p%area(:,i-1,j,k,1)), &
            norm2(p%area(:,i,j,k,2)),norm2(p%area(:,i,j-1,k,2)),norm2(p%area(:,i,j,k-1,3)))
          if(face_max<=0.or.p%legacy_volume(i,j,k)<=0) return
          p%legacy_length(i,j,k)=p%legacy_volume(i,j,k)/face_max
          ! The original revolution grid has equal lower/upper K-face areas.
          ! Cube does not. Preserve legacy_length only as an oracle quantity;
          ! the candidate cube CFL length includes both K faces explicitly.
          face_max=max(face_max,norm2(p%area(:,i,j,k,3)))
          p%cfl_length(i,j,k)=p%legacy_volume(i,j,k)/face_max
        enddo
      enddo
    enddo
    if(minval(p%volume)<=0) return
    status=0
  end subroutine
  pure subroutine face_vertices(p,axis,idx,v,indices)
    class(cube_patch),intent(in)::p
    integer,intent(in)::axis,idx(3)
    real(real64),intent(out)::v(3,4)
    integer,optional,intent(out)::indices(3,4)
    integer::q(3,4),a,b,z
    a=mod(axis,3)+1;b=mod(axis+1,3)+1
    q=spread(idx,dim=2,ncopies=4)
    q(a,1)=idx(a)-1;q(b,1)=idx(b)-1
    q(b,2)=idx(b)-1
    q(a,4)=idx(a)-1
    do z=1,4
      v(:,z)=p%x(:,q(1,z),q(2,z),q(3,z))
    enddo
    if(present(indices)) indices=q
  end subroutine
  pure subroutine edge_ends(p,axis,idx,a,b,indices)
    class(cube_patch),intent(in)::p
    integer,intent(in)::axis,idx(3)
    real(real64),intent(out)::a(3),b(3)
    integer,optional,intent(out)::indices(3,2)
    integer::q(3)
    q=idx;q(axis)=q(axis)-1
    a=p%x(:,q(1),q(2),q(3));b=p%x(:,idx(1),idx(2),idx(3))
    if(present(indices)) then
      indices(:,1)=q;indices(:,2)=idx
    endif
  end subroutine
  pure subroutine patch_entity(p,kind,axis,idx,key,orientation,status)
    class(cube_patch),intent(in)::p
    integer,intent(in)::kind,axis,idx(3) ! 2=edge, 4=face
    integer(int64),intent(out)::key(17)
    integer,intent(out)::orientation,status
    integer::indices(3,4),z,st
    integer(int64)::vertices(4,4)
    real(real64)::v(3,4),a(3),b(3)
    status=1;key=0;orientation=0
    if(axis<1.or.axis>3) return
    if(kind==4) then
      call p%face_vertices(axis,idx,v,indices)
    else if(kind==2) then
      call p%edge_ends(axis,idx,a,b,indices(:,:2))
    else
      return
    endif
    do z=1,kind
      call cube_vertex_key(p%layout%face,p%layout%global, &
        p%layout%first+indices(:,z),vertices(:,z),st)
      if(st/=0) return
    enddo
    call cube_entity_key(vertices(:,:kind),key,orientation);status=0
  end subroutine
  subroutine cube_curl_update(flux,emf,n,dt)
    integer,intent(in)::n(3)
    real(real64),intent(inout)::flux(0:,0:,0:,:)
    real(real64),intent(in)::emf(0:,0:,0:,:),dt
    integer::i,j,k
    do k=1,n(3)
      do j=1,n(2)
        do i=0,n(1)
          flux(i,j,k,1)=flux(i,j,k,1)-dt*(emf(i,j,k,3)-emf(i,j-1,k,3) &
            -emf(i,j,k,2)+emf(i,j,k-1,2))
        enddo
      enddo
    enddo
    do k=1,n(3)
      do j=0,n(2)
        do i=1,n(1)
          flux(i,j,k,2)=flux(i,j,k,2)-dt*(emf(i,j,k,1)-emf(i,j,k-1,1) &
            -emf(i,j,k,3)+emf(i-1,j,k,3))
        enddo
      enddo
    enddo
    do k=0,n(3)
      do j=1,n(2)
        do i=1,n(1)
          flux(i,j,k,3)=flux(i,j,k,3)-dt*(emf(i,j,k,2)-emf(i-1,j,k,2) &
            -emf(i,j,k,1)+emf(i,j-1,k,1))
        enddo
      enddo
    enddo
  end subroutine
  pure function cube_divergence(flux,n) result(div)
    integer,intent(in)::n(3)
    real(real64),intent(in)::flux(0:,0:,0:,:)
    real(real64)::div(n(1),n(2),n(3))
    integer::i,j,k
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          div(i,j,k)=flux(i,j,k,1)-flux(i-1,j,k,1)+flux(i,j,k,2)-flux(i,j-1,k,2) &
            +flux(i,j,k,3)-flux(i,j,k-1,3)
        enddo
      enddo
    enddo
  end function

  pure subroutine cube_inverse3(a,inverse,status)
    real(real64),intent(in)::a(3,3)
    real(real64),intent(out)::inverse(3,3)
    integer,intent(out)::status
    real(real64)::det,scale
    status=1;inverse=0;scale=maxval(abs(a))
    det=dot_product(a(:,1),cube_cross(a(:,2),a(:,3)))
    if(.not.all(ieee_is_finite(a)).or.scale<=tiny(scale)) return
    if(abs(det)<=64*epsilon(det)*scale**3) return
    inverse(1,:)=cube_cross(a(:,2),a(:,3))/det
    inverse(2,:)=cube_cross(a(:,3),a(:,1))/det
    inverse(3,:)=cube_cross(a(:,1),a(:,2))/det
    status=0
  end subroutine
  pure subroutine cube_cell_statistics(v,centre,covariance,volume,status)
    ! Exact degree2 physical polynomial moments of a trilinear hexahedron:
    ! 3-point tensor Gauss integrates degree<=5 per logical coordinate.
    ! Local coordinates avoid subtracting nearly equal planet-scale squares.
    real(real64),intent(in)::v(3,0:1,0:1,0:1)
    real(real64),intent(out)::centre(3),covariance(3,3),volume
    integer,intent(out)::status
    real(real64)::g(3),weight(3),u,w,t,x(3),dr(3),da(3),db(3),p(3),origin(3)
    real(real64)::first(3),second(3,3),jac,q
    integer::i,j,k,a,b,c,d,e
    status=1;centre=0;covariance=0;volume=0;first=0;second=0
    origin=sum(sum(sum(v,dim=4),dim=3),dim=2)/8
    g=[.5_real64-sqrt(3._real64/5)/2,.5_real64,.5_real64+sqrt(3._real64/5)/2]
    weight=[5._real64/18,4._real64/9,5._real64/18]
    do k=1,3
      t=g(k)
      do j=1,3
        w=g(j)
        do i=1,3
          u=g(i);x=0;dr=0;da=0;db=0
          do c=0,1
            do b=0,1
              do a=0,1
                p=v(:,a,b,c)-origin
                x=x+merge(u,1-u,a==1)*merge(w,1-w,b==1)*merge(t,1-t,c==1)*p
                dr=dr+(2*a-1)*merge(w,1-w,b==1)*merge(t,1-t,c==1)*p
                da=da+merge(u,1-u,a==1)*(2*b-1)*merge(t,1-t,c==1)*p
                db=db+merge(u,1-u,a==1)*merge(w,1-w,b==1)*(2*c-1)*p
              enddo
            enddo
          enddo
          jac=dot_product(dr,cube_cross(da,db))
          if(.not.ieee_is_finite(jac).or.jac<=0) return
          q=weight(i)*weight(j)*weight(k)*jac
          volume=volume+q;first=first+q*x
          do e=1,3
            do d=1,3
              second(d,e)=second(d,e)+q*x(d)*x(e)
            enddo
          enddo
        enddo
      enddo
    enddo
    first=first/volume;centre=origin+first;covariance=second/volume
    do e=1,3
      do d=1,3
        covariance(d,e)=covariance(d,e)-first(d)*first(e)
      enddo
    enddo
    status=0
  end subroutine
  pure function cube_polynomial_moments(centre,covariance,origin,transform) result(moments)
    real(real64),intent(in)::centre(3),covariance(3,3),origin(3),transform(3,3)
    real(real64)::moments(10),x(3),second(3,3)
    integer::i,j
    x=matmul(transform,centre-origin)
    second=matmul(transform,matmul(covariance,transpose(transform)))
    do j=1,3
      do i=1,3
        second(i,j)=second(i,j)+x(i)*x(j)
      enddo
    enddo
    moments=[1._real64,x,second(1,1),second(2,2),second(3,3),second(1,2),second(1,3),second(2,3)]
  end function

  subroutine cube_bzz_calibration(p,correction,max_condition,status)
    ! Geometry-only response of original BZZ to three uniform Cartesian Bs.
    ! Cache at setup; multiply the raw BZZ output by correction in the future
    ! cube adapter. Do not change shared CT fluxes or the default LFM operator.
    type(cube_patch),intent(in)::p
    real(real64),allocatable,intent(out)::correction(:,:,:,:,:)
    real(real64),intent(out)::max_condition
    integer,intent(out)::status
    integer::i,j,k,d,c,n(3),idx(3)
    real(real64)::v(3,4),a(3),b(3),delta(3,3),areas(3,3),centre(3),response(3,3),inverse(3,3)
    real(real64)::vol,det,scale
    n=p%layout%local;status=1;max_condition=0
    allocate(correction(3,3,n(1),n(2),n(3)))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          do d=1,3
            idx=[i,j,k];call p%face_vertices(d,idx,v)
            centre=sum(v,dim=2)/4;a=cube_quad(v)
            idx(d)=idx(d)-1;call p%face_vertices(d,idx,v);b=cube_quad(v)
            delta(:,d)=centre-sum(v,dim=2)/4;areas(:,d)=(a+b)/2
          enddo
          vol=abs(dot_product(delta(:,1),cube_cross(delta(:,2),delta(:,3))))
          if(vol<=tiny(vol)) return
          response=0
          do d=1,3
            do c=1,3
              response(:,c)=response(:,c)+delta(:,d)*areas(c,d)/vol
            enddo
          enddo
          scale=maxval(sum(abs(response),dim=2))
          det=dot_product(response(:,1),cube_cross(response(:,2),response(:,3)))
          if(.not.all(ieee_is_finite(response))) return
          if(abs(det)<=64*epsilon(det)*scale**3) return
          inverse(1,:)=cube_cross(response(:,2),response(:,3))/det
          inverse(2,:)=cube_cross(response(:,3),response(:,1))/det
          inverse(3,:)=cube_cross(response(:,1),response(:,2))/det
          max_condition=max(max_condition,scale*maxval(sum(abs(inverse),dim=2)))
          correction(:,:,i,j,k)=inverse
        enddo
      enddo
    enddo
    status=0
  end subroutine
end module
