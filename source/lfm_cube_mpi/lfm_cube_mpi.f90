! Persistent sparse exchange for shared integral faces and EMF edges.
! Plans route integer topology keys once; timesteps exchange values only.
! A canonical owner supplies a single value (no post-CT magnetic averaging).
module lfm_cube_mpi
  use iso_fortran_env,only:real64,int64
  use iso_c_binding,only:c_int,c_double
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use mpi_f08
  use lfm_cube_topology
  use lfm_cube_mesh
  implicit none
  private
  real(real64),allocatable,save::halo_ordered_values(:,:)
  integer(int64),save::halo_ordered_calls=0
  public :: cube_shared_plan, cube_halo_plan, polynomial_weights
  type cube_shared_plan
    type(MPI_Comm) :: comm=MPI_COMM_NULL
    integer :: entries=0,received=0,ranks=0,kind=0,components=1
    integer,allocatable :: index(:,:),orientation(:),permutation(:),leader(:),partner(:)
    integer,allocatable :: sc(:),rc(:),sd(:),rd(:)
    integer,allocatable :: scv(:),rcv(:),sdv(:),rdv(:)
    real(real64),allocatable :: send(:),recv(:),reply(:),back(:)
    ! Derived topology only: actual cell across each conforming patch face.
    ! These references are not point-sampled ghosts or checkpoint state.
    logical,allocatable :: cell_pair_mask(:)
    integer,allocatable :: cell_pair_lookup(:,:,:,:)
  contains
    procedure :: init => shared_init
    procedure :: sync => shared_sync
    procedure :: sync_fields => shared_sync_fields
    procedure :: prepare_cell_pairs => shared_prepare_cell_pairs
    procedure :: exchange_base => shared_exchange_base
    procedure :: close => shared_close
  end type
  type cube_halo_plan
    ! One setup of unique donor cells; value-only exchanges thereafter.
    ! States are Cartesian cell averages, not face-normal CT fluxes.
    type(MPI_Comm) :: comm=MPI_COMM_NULL
    type(cube_layout) :: layout
    integer :: ng=0,nvars=0,targets=0,donors=0,received=0
    integer,allocatable :: index(:,:),reference(:,:),terms(:),permutation(:),local_donor(:,:)
    integer,allocatable :: sc(:),rc(:),sd(:),rd(:),scv(:),rcv(:),sdv(:),rdv(:)
    real(real64),allocatable :: weight(:,:),send(:),recv(:),samples(:,:)
    real(real64) :: maximum_moment_residual=0,maximum_weight_norm=0
    real(real64) :: profile_seconds(4)=0 ! pack, MPI, sample unpack, interpolation
  contains
    procedure :: init => halo_init
    procedure :: sync => halo_sync
    procedure :: close => halo_close
  end type
  interface
    integer(c_int) function halo_ordered8(nt,nd,ldr,ldw,terms,reference,weight,samples,positive,result,limited) bind(C,name='cube_yy_halo_ordered8')
      import c_int,c_double
      integer(c_int),value::nt,nd,ldr,ldw
      integer(c_int),intent(in)::terms(*),reference(*),positive(*)
      real(c_double),intent(in)::weight(*),samples(*)
      real(c_double),intent(out)::result(*)
      integer(c_int),intent(out)::limited
    end function
  end interface
contains
  subroutine halo_init(p,layout,meridian,ng,nvars,comm,status)
    class(cube_halo_plan),intent(out)::p
    type(cube_layout),intent(in)::layout
    type(cube_meridian),intent(in)::meridian
    integer,intent(in)::ng,nvars
    type(MPI_Comm),intent(in)::comm
    integer,intent(out)::status
    integer,allocatable::cells(:,:),owners(:),cursor(:),packed(:,:),incoming(:,:)
    integer::st,bad,i,q,pos,ranks,local(3),owner
    call build_halo_local(p,layout,meridian,ng,nvars,cells,status)
    call MPI_Allreduce(status,bad,1,MPI_INTEGER,MPI_MAX,comm,st)
    status=1
    if(st/=MPI_SUCCESS.or.bad/=0) return
    call MPI_Comm_dup(comm,p%comm,st)
    if(st/=MPI_SUCCESS) return
    call MPI_Comm_size(p%comm,ranks,st)
    allocate(owners(p%donors),p%sc(ranks),p%rc(ranks),p%sd(ranks),p%rd(ranks),cursor(ranks))
    p%sc=0
    do i=1,p%donors
      call cube_owner(cells(1,i),cells(2:,i)-1,layout%global,layout%dims,owner,local,st)
      if(st/=0) error stop 'halo owner invariant'
      owners(i)=owner+1;p%sc(owner+1)=p%sc(owner+1)+1
    enddo
    call MPI_Alltoall(p%sc,1,MPI_INTEGER,p%rc,1,MPI_INTEGER,p%comm,st)
    if(st/=MPI_SUCCESS) return
    p%sd(1)=0;p%rd(1)=0
    do q=2,ranks
      p%sd(q)=p%sd(q-1)+p%sc(q-1);p%rd(q)=p%rd(q-1)+p%rc(q-1)
    enddo
    p%received=sum(p%rc)
    allocate(p%permutation(p%donors),packed(4,p%donors),incoming(4,p%received))
    cursor=p%sd
    do i=1,p%donors
      q=owners(i);pos=cursor(q)+1;cursor(q)=pos
      packed(:,pos)=cells(:,i);p%permutation(pos)=i
    enddo
    call MPI_Alltoallv(packed,4*p%sc,4*p%sd,MPI_INTEGER,incoming,4*p%rc,4*p%rd,MPI_INTEGER,p%comm,st)
    if(st/=MPI_SUCCESS) return
    allocate(p%local_donor(3,p%received));bad=0
    do i=1,p%received
      local=incoming(2:,i)-layout%first
      if(incoming(1,i)/=layout%face.or.any(local<1).or.any(local>layout%local)) bad=1
      p%local_donor(:,i)=local
    enddo
    call MPI_Allreduce(bad,q,1,MPI_INTEGER,MPI_MAX,p%comm,st)
    if(st/=MPI_SUCCESS.or.q/=0) return
    p%scv=nvars*p%sc;p%rcv=nvars*p%rc;p%sdv=nvars*p%sd;p%rdv=nvars*p%rd
    allocate(p%send(nvars*p%received),p%recv(nvars*p%donors),p%samples(nvars,p%donors))
    status=0
  end subroutine

  subroutine build_halo_local(p,layout,m,ng,nvars,cells,status)
    class(cube_halo_plan),intent(inout)::p
    type(cube_layout),intent(in)::layout
    type(cube_meridian),intent(in)::m
    integer,intent(in)::ng,nvars
    integer,allocatable,intent(out)::cells(:,:)
    integer,intent(out)::status
    integer(int64),allocatable::keys(:,:)
    integer,allocatable::order(:),inverse(:),req(:,:),home(:,:)
    real(real64),allocatable::centres(:,:),covs(:,:,:)
    real(real64)::v(3,0:1,0:1,0:1),centre(3),cov(3,3),vol,a,b,origin(3),transform(3,3),delta(3,3)
    real(real64)::matrix(27,10),target(10),weights(27),residual,normw,pi,radial_t,direction(3)
    integer::n(3),g(3),idx(3),start(3),f,i,j,k,h,t,z,cap,requests,d,a0,b0,c0,st,unique
    status=1;n=layout%local;pi=acos(-1._real64)
    if(ng<1.or.nvars<1.or.any(layout%global<3).or.layout%global(1)/=m%ni) return
    if(.not.allocated(m%xy).or.any(n<=0).or.2*ng>=layout%global(2)) return
    p%layout=layout;p%ng=ng;p%nvars=nvars
    cap=(n(1)+2*ng)*(n(2)+2*ng)*(n(3)+2*ng)-product(n)
    allocate(p%index(3,cap),p%terms(cap),p%reference(27,cap),p%weight(27,cap), &
      req(4,27*cap),home(4,cap))
    requests=0;t=0;p%reference=0;p%weight=0
    do k=1-ng,n(3)+ng
      do j=1-ng,n(2)+ng
        do i=1-ng,n(1)+ng
          if(i>=1.and.i<=n(1).and.j>=1.and.j<=n(2).and.k>=1.and.k<=n(3)) cycle
          g=layout%first+[i,j,k]
          ! Radial physical boundary guards belong to the physical BC, not MPI.
          if(g(1)<1.or.g(1)>layout%global(1)) cycle
          t=t+1;p%index(:,t)=[i,j,k]
          if(all(g(2:)>=1).and.all(g(2:)<=layout%global(2:))) then
            requests=requests+1;req(:,requests)=[layout%face,g]
            p%terms(t)=1;p%reference(1,t)=requests;p%weight(1,t)=1
            cycle
          endif
          call m%cell_vertices(layout%face,layout%global(2),g,v,st)
          if(st/=0) return
          call cube_cell_statistics(v,centre,cov,vol,st)
          if(st/=0) return
          if(m%shell%ns>0) then
            call m%reference(centre,radial_t,direction,st)
            if(st/=0) return
            call cube_inverse(direction,f,a,b,st)
          else
            call cube_inverse(centre,f,a,b,st)
          endif
          if(st/=0) return
          idx=[g(1),int(floor((a/(pi/4)+1)*layout%global(2)/2))+1, &
            int(floor((b/(pi/4)+1)*layout%global(3)/2))+1]
          if(m%shell%ns>0) idx(1)=int(floor(radial_t*layout%global(1)))+1
          idx=max(1,min(layout%global,idx));home(:,t)=[f,idx]
          start=max(1,min(layout%global-2,idx-1));z=0
          do c0=0,2
            do b0=0,2
              do a0=0,2
                requests=requests+1;z=z+1;req(:,requests)=[f,start+[a0,b0,c0]]
                p%reference(z,t)=requests
              enddo
            enddo
          enddo
          p%terms(t)=27
        enddo
      enddo
    enddo
    p%targets=t
    allocate(keys(4,requests),order(requests),inverse(requests))
    keys=int(req(:,:requests),int64);order=[(i,i=1,requests)]
    call sort_keys(keys,order,1,requests)
    unique=0
    do i=1,requests
      if(i==1) then
        unique=unique+1
      else if(any(keys(:,order(i))/=keys(:,order(i-1)))) then
        unique=unique+1
      endif
      inverse(order(i))=unique
    enddo
    p%donors=unique;allocate(cells(4,unique));unique=0
    do i=1,requests
      if(inverse(order(i))>unique) then
        unique=unique+1;cells(:,unique)=req(:,order(i))
      endif
    enddo
    deallocate(req,keys,order)
    allocate(centres(3,unique),covs(3,3,unique))
    do i=1,unique
      call m%cell_vertices(cells(1,i),layout%global(2),cells(2:,i),v,st)
      if(st/=0) return
      call cube_cell_statistics(v,centres(:,i),covs(:,:,i),vol,st)
      if(st/=0) return
    enddo
    do h=1,t
      d=p%terms(h);p%reference(:d,h)=inverse(p%reference(:d,h))
      if(d==1) cycle
      call m%cell_vertices(home(1,h),layout%global(2),home(2:,h),v,st)
      if(st/=0) return
      delta(:,1)=(sum(sum(v(:,1,:,:),dim=3),dim=2)-sum(sum(v(:,0,:,:),dim=3),dim=2))/4
      delta(:,2)=(sum(sum(v(:,:,1,:),dim=3),dim=2)-sum(sum(v(:,:,0,:),dim=3),dim=2))/4
      delta(:,3)=(sum(sum(v(:,:,:,1),dim=3),dim=2)-sum(sum(v(:,:,:,0),dim=3),dim=2))/4
      call cube_inverse3(delta,transform,st)
      if(st/=0) return
      origin=sum(sum(sum(v,dim=4),dim=3),dim=2)/8
      do z=1,27
        i=p%reference(z,h)
        matrix(z,:)=cube_polynomial_moments(centres(:,i),covs(:,:,i),origin,transform)
      enddo
      call m%cell_vertices(layout%face,layout%global(2),layout%first+p%index(:,h),v,st)
      if(st/=0) return
      call cube_cell_statistics(v,centre,cov,vol,st)
      if(st/=0) return
      target=cube_polynomial_moments(centre,cov,origin,transform)
      call polynomial_weights(matrix,target,weights,residual,st)
      if(st/=0) return
      p%weight(:,h)=weights;normw=sum(abs(weights))
      p%maximum_moment_residual=max(p%maximum_moment_residual,residual)
      p%maximum_weight_norm=max(p%maximum_weight_norm,normw)
    enddo
    status=0
  end subroutine

  subroutine polynomial_weights(matrix,target,weights,residual,status)
    ! Two-pass modified Gram-Schmidt QR; no squared-condition normal equations.
    real(real64),intent(in)::matrix(27,10),target(10)
    real(real64),intent(out)::weights(27),residual
    integer,intent(out)::status
    real(real64)::q(27,10),r(10,10),v(27),z(10),c,scale
    integer::i,j,pass
    status=1;r=0;q=0;scale=maxval(abs(matrix))
    do j=1,10
      v=matrix(:,j)
      do pass=1,2
        do i=1,j-1
          c=dot_product(q(:,i),v);r(i,j)=r(i,j)+c;v=v-c*q(:,i)
        enddo
      enddo
      r(j,j)=sqrt(sum(v*v))
      if(r(j,j)<=256*epsilon(scale)*scale) return
      q(:,j)=v/r(j,j)
    enddo
    z=0
    do j=1,10
      z(j)=(target(j)-dot_product(r(:j-1,j),z(:j-1)))/r(j,j)
    enddo
    weights=matmul(q,z)
    residual=maxval(abs(matmul(transpose(matrix),weights)-target))/max(1._real64,maxval(abs(target)))
    if(residual>1.e-10_real64) return
    status=0
  end subroutine

  subroutine halo_sync(p,values,status,bounded,positive,limited)
    class(cube_halo_plan),intent(inout)::p
    real(real64),intent(inout)::values(1-p%ng:,1-p%ng:,1-p%ng:,:)
    integer,intent(out)::status
    logical,optional,intent(in)::bounded(:),positive(:)
    integer,optional,intent(out)::limited
    integer::i,j,k,h,z,st,idx(3),d,count
    logical::need_lo(p%nvars),need_hi(p%nvars)
    real(real64)::value,lo,hi,tick
    logical::use_ordered,oracle,words_ok
    integer(c_int)::ordered_status,ordered_limited,posmask(8)
    integer::cache_status,unit,ios,rank
    character(128)::packet_name
    status=1;count=0;tick=MPI_Wtime()
    if(p%comm==MPI_COMM_NULL.or.size(values,4)/=p%nvars) return
    if(present(bounded)) then
      if(size(bounded)/=p%nvars) return
    endif
    if(present(positive)) then
      if(size(positive)/=p%nvars) return
    endif
    ! The caller mask is constant over all targets. Same-face targets
    ! have exactly one donor with unit weight: copy without the generic
    ! 27-point gather/reduction machinery. Cross-face DOT_PRODUCT order and
    ! all donor positivity/bounding rules remain unchanged.
    need_lo=.false.;need_hi=.false.
    if(present(bounded)) then
      need_lo=need_lo.or.bounded;need_hi=bounded
    endif
    do i=1,p%received
      idx=p%local_donor(:,i)
      do z=1,p%nvars
        p%send(p%nvars*(i-1)+z)=values(idx(1),idx(2),idx(3),z)
      enddo
    enddo
    p%profile_seconds(1)=p%profile_seconds(1)+MPI_Wtime()-tick;tick=MPI_Wtime()
    call sparse_values(p%comm,p%scv,p%rcv,p%sdv,p%rdv,p%send,p%recv,.true.,st)
    p%profile_seconds(2)=p%profile_seconds(2)+MPI_Wtime()-tick;tick=MPI_Wtime()
    if(st/=0) return
    do i=1,p%donors
      j=p%permutation(i)
      p%samples(:,j)=p%recv(p%nvars*(i-1)+1:p%nvars*i)
    enddo
    p%profile_seconds(3)=p%profile_seconds(3)+MPI_Wtime()-tick;tick=MPI_Wtime()
    if(.not.all(ieee_is_finite(p%samples))) return
    ! Unique donors are exactly the union of all target stencils. Check
    ! positive channels once; retain the original stencil MIN only where
    ! bounding or the original negative-extrapolation fallback needs it.
    if(present(positive)) then
      do z=1,p%nvars
        if(positive(z)) then
          if(any(p%samples(z,:)<=0)) return ! Invalid donors fail; no floor.
        endif
      enddo
    endif
    use_ordered=.false.;oracle=.false.
    if(p%nvars==8.and..not.any(need_lo).and..not.any(need_hi))then
      cache_status=0
      if(allocated(halo_ordered_values))then
        if(any(shape(halo_ordered_values)/=[8,p%targets]))deallocate(halo_ordered_values)
      endif
      if(.not.allocated(halo_ordered_values))allocate(halo_ordered_values(8,p%targets),stat=cache_status)
      if(cache_status==0)then
        posmask=0
        if(present(positive))posmask=merge(1_c_int,0_c_int,positive)
        ordered_status=halo_ordered8(int(p%targets,c_int),int(p%donors,c_int),int(size(p%reference,1),c_int), &
          int(size(p%weight,1),c_int),p%terms,p%reference,p%weight,p%samples,posmask,halo_ordered_values,ordered_limited)
        use_ordered=ordered_status==0
        if(use_ordered)then
          halo_ordered_calls=halo_ordered_calls+1
          ! Initial cases and adjacent power-of-two calls cover current and
          ! predictor masks during real evolution and own restart. Keep the
          ! original Fortran DOT_PRODUCT as a word oracle on those calls.
          oracle=halo_ordered_calls<=8
          if(halo_ordered_calls>1)oracle=oracle.or.iand(halo_ordered_calls,halo_ordered_calls-1)==0.or. &
            iand(halo_ordered_calls-1,halo_ordered_calls-2)==0
        endif
      endif
    endif
    if(.not.use_ordered.or.oracle)then
    do h=1,p%targets
      idx=p%index(:,h);d=p%terms(h)
      if(d==1) then
        do z=1,p%nvars
          value=p%samples(z,p%reference(1,h))
          if(present(positive)) then
            if(positive(z).and.value<=0) return
          endif
          values(idx(1),idx(2),idx(3),z)=value
        enddo
        cycle
      endif
      do z=1,p%nvars
        if(need_lo(z)) lo=minval(p%samples(z,p%reference(:d,h)))
        if(need_hi(z)) hi=maxval(p%samples(z,p%reference(:d,h)))
        value=dot_product(p%weight(:d,h),p%samples(z,p%reference(:d,h)))
        if(present(positive)) then
          if(positive(z)) then
            if(value<=0) then
              if(.not.need_lo(z)) lo=minval(p%samples(z,p%reference(:d,h)))
              value=lo;count=count+1
            endif
          endif
        endif
        if(present(bounded)) then
          if(bounded(z)) then
            if(value<lo.or.value>hi) count=count+1
            value=max(lo,min(hi,value))
          endif
        endif
        values(idx(1),idx(2),idx(3),z)=value
      enddo
    enddo
    else
      count=int(ordered_limited)
      do h=1,p%targets
        idx=p%index(:,h)
        do z=1,8
          values(idx(1),idx(2),idx(3),z)=halo_ordered_values(z,h)
        enddo
      enddo
    endif
    if(use_ordered.and.oracle)then
      words_ok=count==int(ordered_limited)
      do h=1,p%targets
        idx=p%index(:,h)
        do z=1,8
          if(transfer(values(idx(1),idx(2),idx(3),z),0_int64)/=transfer(halo_ordered_values(z,h),0_int64))words_ok=.false.
        enddo
      enddo
      if(.not.words_ok)then
        call MPI_Comm_rank(p%comm,rank,st)
        write(packet_name,'("halo-ordered-mismatch-r",i3.3,"-call",i0,".bin")')rank,halo_ordered_calls
        open(newunit=unit,file=trim(packet_name),access='stream',form='unformatted',status='new',iostat=ios)
        if(ios==0)then
          write(unit)'YYHALO24',p%targets,p%donors,size(p%reference,1),size(p%weight,1),count,ordered_limited, &
            p%terms,p%reference,p%weight,p%samples,posmask
          do h=1,p%targets
            idx=p%index(:,h);write(unit)values(idx(1),idx(2),idx(3),:)
          enddo
          write(unit)halo_ordered_values
          close(unit)
        endif
        status=1;return
      endif
      call MPI_Comm_rank(p%comm,rank,st)
      write(*,*)'ORDERED_HALO_ACTUAL_WORD_PASS rank/call/targets=',rank,halo_ordered_calls,p%targets
    endif
    p%profile_seconds(4)=p%profile_seconds(4)+MPI_Wtime()-tick
    if(present(limited)) limited=count
    status=0
  end subroutine
  subroutine halo_close(p)
    class(cube_halo_plan),intent(inout)::p
    integer::st
    if(p%comm/=MPI_COMM_NULL) call MPI_Comm_free(p%comm,st)
  end subroutine

  subroutine shared_init(p,mesh,kind,comm,status,components,reconcile)
    class(cube_shared_plan),intent(out)::p
    type(cube_patch),intent(in)::mesh
    integer,intent(in)::kind
    type(MPI_Comm),intent(in)::comm
    integer,intent(out)::status
    integer,optional,intent(in)::components
    logical,optional,intent(in)::reconcile
    integer::n(3),cap,i,j,k,d,a,b,m,st,q,z,bucket,pos,first,last,chosen,other,perface
    integer,allocatable::owners(:),cursor(:),order(:),kc(:),kr(:),kd(:),kdr(:),source_face(:)
    integer(int64),allocatable::keys(:,:),packed(:,:),received(:,:)
    integer(int64)::key(17),h
    logical::include,pair_average
    status=1;n=mesh%layout%local
    pair_average=.false.
    if(present(reconcile)) pair_average=reconcile
    if(pair_average.and.kind/=4) return
    if(present(components)) p%components=components
    if(p%components<1) return
    if(kind/=2.and.kind/=4) return
    if(any(n<=0)) return
    call MPI_Comm_dup(comm,p%comm,st)
    if(st/=MPI_SUCCESS) return
    call MPI_Comm_size(p%comm,p%ranks,st)
    if(st/=MPI_SUCCESS) return
    p%kind=kind
    if(kind==4) then
      cap=2*(n(1)*n(2)+n(2)*n(3)+n(1)*n(3))
    else
      cap=4*(n(1)*n(2)+n(2)*n(3)+n(1)*n(3))
    endif
    allocate(p%index(4,cap),p%orientation(cap),keys(17,cap),owners(cap))
    m=0
    do d=1,3
      a=mod(d,3)+1;b=mod(d+1,3)+1
      do k=0,n(3)
        do j=0,n(2)
          do i=0,n(1)
            block
              integer::v(3)
              v=[i,j,k]
              if(kind==4) then
                include=(v(d)==0.or.v(d)==n(d)).and.v(a)>0.and.v(b)>0
              else
                include=v(d)>0.and.(v(a)==0.or.v(a)==n(a).or.v(b)==0.or.v(b)==n(b))
              endif
            end block
            if(.not.include) cycle
            m=m+1
            if(m>cap) error stop 'cube shared capacity'
            call mesh%entity(kind,d,[i,j,k],key,q,st)
            if(st/=0) return
            p%index(:,m)=[i,j,k,d];p%orientation(m)=q;keys(:,m)=key
            h=0
            do z=1,17
              h=modulo(131_int64*h+key(z),int(p%ranks,int64))
            enddo
            owners(m)=int(h)+1
          enddo
        enddo
      enddo
    enddo
    p%entries=m
    allocate(p%sc(p%ranks),p%rc(p%ranks),p%sd(p%ranks),p%rd(p%ranks),cursor(p%ranks))
    p%sc=0
    do i=1,m
      p%sc(owners(i))=p%sc(owners(i))+1
    enddo
    call MPI_Alltoall(p%sc,1,MPI_INTEGER,p%rc,1,MPI_INTEGER,p%comm,st)
    if(st/=MPI_SUCCESS) return
    p%sd(1)=0;p%rd(1)=0
    do i=2,p%ranks
      p%sd(i)=p%sd(i-1)+p%sc(i-1);p%rd(i)=p%rd(i-1)+p%rc(i-1)
    enddo
    p%received=sum(p%rc)
    allocate(p%permutation(m),packed(17,m),received(17,p%received))
    cursor=p%sd
    do i=1,m
      bucket=owners(i);pos=cursor(bucket)+1;cursor(bucket)=pos
      packed(:,pos)=keys(:,i);p%permutation(pos)=i
    enddo
    kc=17*p%sc;kr=17*p%rc;kd=17*p%sd;kdr=17*p%rd
    call MPI_Alltoallv(packed,kc,kd,MPI_INTEGER8,received,kr,kdr,MPI_INTEGER8,p%comm,st)
    if(st/=MPI_SUCCESS) return
    allocate(order(p%received),p%leader(p%received),p%partner(p%received),source_face(p%received))
    perface=product(mesh%layout%dims)
    do i=1,p%ranks
      source_face(p%rd(i)+1:p%rd(i)+p%rc(i))=(i-1)/perface
    enddo
    order=[(i,i=1,p%received)]
    call sort_keys(received,order,1,p%received)
    first=1
    do while(first<=p%received)
      last=first
      do while(last<p%received)
        if(any(received(:,order(last+1))/=received(:,order(first)))) exit
        last=last+1
      enddo
      ! Received segments are in source-rank order. A minimum index chooses
      ! a deterministic source without summing or changing matching values.
      chosen=minval(order(first:last))
      other=chosen
      if(pair_average.and.last-first==1) then
        j=maxval(order(first:last))
        ! GAMERA-Cube's MPI fluid seam policy: reconcile the two oriented
        ! reconstructions. Internal same-patch faces keep exact ownership.
        if(source_face(chosen)/=source_face(j)) other=j
      endif
      do i=first,last
        p%leader(order(i))=chosen
        p%partner(order(i))=other
      enddo
      if(kind==4.and.last-first+1>2) error stop 'non-manifold cube face'
      if(kind==2.and.last-first+1>4) error stop 'non-manifold cube edge'
      first=last+1
    enddo
    p%scv=p%components*p%sc;p%rcv=p%components*p%rc
    p%sdv=p%components*p%sd;p%rdv=p%components*p%rd
    allocate(p%send(p%components*m),p%back(p%components*m), &
      p%recv(p%components*p%received),p%reply(p%components*p%received))
    status=0
  end subroutine
  recursive subroutine sort_keys(keys,order,lo,hi)
    integer(int64),intent(in)::keys(:,:)
    integer,intent(inout)::order(:)
    integer,intent(in)::lo,hi
    integer::i,j,t
    integer(int64)::pivot(size(keys,1))
    if(lo>=hi) return
    i=lo;j=hi;pivot=keys(:,order(lo+(hi-lo)/2))
    do
      do while(cube_key_less(keys(:,order(i)),pivot))
        i=i+1
      enddo
      do while(cube_key_less(pivot,keys(:,order(j))))
        j=j-1
      enddo
      if(i>j) exit
      t=order(i);order(i)=order(j);order(j)=t;i=i+1;j=j-1
      if(i>j) exit
    enddo
    if(lo<j) call sort_keys(keys,order,lo,j)
    if(i<hi) call sort_keys(keys,order,i,hi)
  end subroutine
  subroutine sparse_values(comm,sc,rc,sd,rd,send,recv,reverse,status)
    type(MPI_Comm),intent(in)::comm
    integer,intent(in)::sc(:),rc(:),sd(:),rd(:)
    real(real64),contiguous,intent(in),asynchronous::send(:)
    real(real64),contiguous,intent(out),asynchronous::recv(:)
    logical,intent(in)::reverse
    integer,intent(out)::status
    type(MPI_Request)::requests(2*size(sc))
    integer::q,count,st,tag
    count=0;tag=111
    if(reverse) tag=112
    do q=1,size(sc)
      if(.not.reverse) then
        if(rc(q)==0) cycle
        count=count+1
        call MPI_Irecv(recv(rd(q)+1:),rc(q),MPI_DOUBLE_PRECISION,q-1,tag,comm,requests(count),st)
      else
        if(sc(q)==0) cycle
        count=count+1
        call MPI_Irecv(recv(sd(q)+1:),sc(q),MPI_DOUBLE_PRECISION,q-1,tag,comm,requests(count),st)
      endif
      if(st/=MPI_SUCCESS) error stop 'cube MPI_Irecv'
    enddo
    do q=1,size(sc)
      if(.not.reverse) then
        if(sc(q)==0) cycle
        count=count+1
        call MPI_Isend(send(sd(q)+1:),sc(q),MPI_DOUBLE_PRECISION,q-1,tag,comm,requests(count),st)
      else
        if(rc(q)==0) cycle
        count=count+1
        call MPI_Isend(send(rd(q)+1:),rc(q),MPI_DOUBLE_PRECISION,q-1,tag,comm,requests(count),st)
      endif
      if(st/=MPI_SUCCESS) error stop 'cube MPI_Isend'
    enddo
    call MPI_Waitall(count,requests,MPI_STATUSES_IGNORE,st)
    status=merge(0,1,st==MPI_SUCCESS)
  end subroutine
  subroutine shared_sync(p,values,status,max_jump)
    class(cube_shared_plan),intent(inout)::p
    real(real64),intent(inout)::values(0:,0:,0:,:)
    integer,intent(out)::status
    real(real64),optional,intent(out)::max_jump
    integer::i,j,q(4),st
    real(real64)::jump
    status=1
    if(p%comm==MPI_COMM_NULL) return
    if(p%components/=1) return
    do i=1,p%entries
      j=p%permutation(i);q=p%index(:,j)
      p%send(i)=p%orientation(j)*values(q(1),q(2),q(3),q(4))
    enddo
    call sparse_values(p%comm,p%sc,p%rc,p%sd,p%rd,p%send,p%recv,.false.,st)
    if(st/=0) return
    do i=1,p%received
      p%reply(i)=p%recv(p%leader(i))
    enddo
    if(present(max_jump)) then
      jump=0
      if(p%received>0) jump=maxval(abs(p%reply-p%recv))
      call MPI_Allreduce(jump,max_jump,1,MPI_DOUBLE_PRECISION,MPI_MAX,p%comm,st)
      if(st/=MPI_SUCCESS) return
    endif
    call sparse_values(p%comm,p%sc,p%rc,p%sd,p%rd,p%reply,p%back,.true.,st)
    if(st/=0) return
    do i=1,p%entries
      j=p%permutation(i);q=p%index(:,j)
      values(q(1),q(2),q(3),q(4))=p%orientation(j)*p%back(i)
    enddo
    status=0
  end subroutine
  subroutine shared_sync_fields(p,values,status,max_jump)
    ! All Cartesian flux channels share one orientation. Opt-in fluid seams
    ! use GAMERA-Cube's two-sided average; CT/EMF and same-patch ownership
    ! retain their original path. No per-step metadata or allocations.
    class(cube_shared_plan),intent(inout)::p
    real(real64),intent(inout)::values(0:,0:,0:,:,:)
    integer,intent(out)::status
    real(real64),optional,intent(out)::max_jump(:)
    integer::i,j,q(4),st,c,offset,source,other
    real(real64)::jump(p%components)
    status=1
    if(p%comm==MPI_COMM_NULL) return
    if(size(values,4)/=3.or.size(values,5)/=p%components) return
    if(present(max_jump)) then
      if(size(max_jump)/=p%components) return
    endif
    do i=1,p%entries
      j=p%permutation(i);q=p%index(:,j);offset=p%components*(i-1)
      do c=1,p%components
        p%send(offset+c)=p%orientation(j)*values(q(1),q(2),q(3),q(4),c)
      enddo
    enddo
    call sparse_values(p%comm,p%scv,p%rcv,p%sdv,p%rdv,p%send,p%recv,.false.,st)
    if(st/=0) return
    jump=0
    do i=1,p%received
      offset=p%components*(i-1);source=p%components*(p%leader(i)-1)
      other=p%components*(p%partner(i)-1)
      do c=1,p%components
        p%reply(offset+c)=p%recv(source+c)
        if(other/=source) p%reply(offset+c)=.5_real64*(p%recv(source+c)+p%recv(other+c))
        if(present(max_jump)) jump(c)=max(jump(c),abs(p%reply(offset+c)-p%recv(offset+c)))
      enddo
    enddo
    if(present(max_jump)) then
      call MPI_Allreduce(jump,max_jump,p%components,MPI_DOUBLE_PRECISION,MPI_MAX,p%comm,st)
      if(st/=MPI_SUCCESS) return
    endif
    call sparse_values(p%comm,p%scv,p%rcv,p%sdv,p%rdv,p%reply,p%back,.true.,st)
    if(st/=0) return
    do i=1,p%entries
      j=p%permutation(i);q=p%index(:,j);offset=p%components*(i-1)
      do c=1,p%components
        values(q(1),q(2),q(3),q(4),c)=p%orientation(j)*p%back(offset+c)
      enddo
    enddo
    status=0
  end subroutine

  subroutine shared_prepare_cell_pairs(p,n,status)
    ! Reuse the exact key/owner plan. Only cross-patch two-member groups
    ! receive a physical neighbour; same-patch and physical boundaries keep
    ! their old point-ghost/base references unchanged.
    class(cube_shared_plan),intent(inout)::p
    integer,intent(in)::n(3)
    integer,intent(out)::status
    integer::i,j,q(4),offset,st
    status=1
    if(p%kind/=4.or.p%components/=8.or.p%comm==MPI_COMM_NULL) return
    if(any(n<=0).or.allocated(p%cell_pair_mask)) return
    allocate(p%cell_pair_mask(p%entries),p%cell_pair_lookup(0:n(1),0:n(2),0:n(3),3))
    p%cell_pair_mask=.false.;p%cell_pair_lookup=0;p%reply=0
    do i=1,p%received
      offset=8*(i-1)
      if(p%leader(i)/=p%partner(i)) p%reply(offset+1)=1
    enddo
    call sparse_values(p%comm,p%scv,p%rcv,p%sdv,p%rdv,p%reply,p%back,.true.,st)
    if(st/=MPI_SUCCESS) return
    do i=1,p%entries
      j=p%permutation(i);q=p%index(:,j);offset=8*(i-1)
      if(q(4)<1.or.q(4)>3) return
      if(any(q(1:3)<0).or.any(q(1:3)>n)) return
      if(p%cell_pair_lookup(q(1),q(2),q(3),q(4))/=0) return
      if(p%back(offset+1)/=0.and.p%back(offset+1)/=1) return
      p%cell_pair_mask(j)=p%back(offset+1)==1
      p%cell_pair_lookup(q(1),q(2),q(3),q(4))=j
    enddo
    status=0
  end subroutine

  subroutine shared_exchange_base(p,g,q,reference,status)
    ! Cartesian primitive states have no flux orientation sign. Exchange
    ! only actual active cells; never rotate, average or clip their values.
    class(cube_shared_plan),intent(inout)::p
    integer,intent(in)::g
    real(real64),intent(in)::q(1-g:,1-g:,1-g:,:)
    real(real64),intent(out)::reference(:,:)
    integer,intent(out)::status
    integer::i,j,idx(4),cell(3),offset,other,st,n(3)
    status=1;n=shape(q(:,:,:,1))-2*g
    if(p%kind/=4.or.p%components/=8.or.size(q,4)/=8) return
    if(.not.allocated(p%cell_pair_mask)) return
    if(any(shape(reference)/=[8,p%entries])) return
    if(any(shape(p%cell_pair_lookup)/=[n+1,3])) return
    do i=1,p%entries
      j=p%permutation(i);idx=p%index(:,j);cell=idx(1:3);offset=8*(i-1)
      if(cell(idx(4))==0) then
        cell(idx(4))=1
      else if(cell(idx(4))/=n(idx(4))) then
        return
      endif
      p%send(offset+1:offset+8)=q(cell(1),cell(2),cell(3),:)
    enddo
    call sparse_values(p%comm,p%scv,p%rcv,p%sdv,p%rdv,p%send,p%recv,.false.,st)
    if(st/=MPI_SUCCESS) return
    do i=1,p%received
      offset=8*(i-1);other=i
      if(p%leader(i)/=p%partner(i)) then
        if(i==p%leader(i)) then
          other=p%partner(i)
        else if(i==p%partner(i)) then
          other=p%leader(i)
        else
          return
        endif
      endif
      p%reply(offset+1:offset+8)=p%recv(8*(other-1)+1:8*other)
    enddo
    call sparse_values(p%comm,p%scv,p%rcv,p%sdv,p%rdv,p%reply,p%back,.true.,st)
    if(st/=MPI_SUCCESS) return
    do i=1,p%entries
      j=p%permutation(i);offset=8*(i-1)
      reference(:,j)=p%back(offset+1:offset+8)
    enddo
    if(.not.all(ieee_is_finite(reference))) return
    if(minval(reference(1,:))<=0.or.minval(reference(5,:))<=0) return
    status=0
  end subroutine

  subroutine shared_close(p)
    class(cube_shared_plan),intent(inout)::p
    integer::st
    if(p%comm/=MPI_COMM_NULL) call MPI_Comm_free(p%comm,st)
  end subroutine
end module
