! Native cap storage/units adapter for the released YY ghost-wall routines.
module lfm_cube_yy_boundary
 use iso_c_binding,only:c_int,c_double,c_int64_t
 use iso_fortran_env,only:int64
 use mpi_f08
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 use lfm_rotating_mi_map,only:mi_cap_weights
 use lfm_cube_yy_native_cap_frame,only:yy_native_cap_frame
 use lfm_cube_exact_cap_lookup,only:exact_cap_configure,exact_cap_weights
 implicit none
 private
 public::yy_boundary_configure,yy_boundary_install,yy_boundary_potential,yy_boundary_ghost
 public::yy_boundary_set_time
 public::yy_boundary_wall_cache_metrics,yy_boundary_volume_domain
 ! Exact memoization of the held-epoch point potential. Current-time B and
 ! corotation in the released C wall are still evaluated on every call.
 public::yy_boundary_cache_enable,yy_boundary_cache_metrics,yy_boundary_collective_refresh
 integer(int64),allocatable,save::point_key(:,:)
 real(c_double),allocatable,save::point_value(:)
 logical,allocatable,save::point_used(:),point_valid(:)
 integer,save::point_count=0
 integer(int64),save::cache_hits=0,cache_misses=0
 logical,save::cache_enabled=.true.
 real(c_double),allocatable,save::cap_x(:,:),cap_y(:,:),north(:,:),south(:,:)
 real(c_double),save::length_m=0,held_time=-1,query_time=-1,spin_omega=0
 logical,save::installed=.false.,exact_cap_usable=.false.
 interface
  integer(c_int) function cube_yy_volume_cache_configure(low,high,rion) bind(C)
   import c_double,c_int
   real(c_double),value::low,high,rion
  end function
  integer(c_int) function cube_yy_volume_cache_refresh(comm,time) bind(C)
   import c_double,c_int
   integer(c_int),value::comm
   real(c_double),value::time
  end function
  integer(c_int) function cube_yy_volume_potential(point,value) bind(C)
   import c_double,c_int
   real(c_double),intent(in)::point(3)
   real(c_double),intent(out)::value
  end function
  subroutine cube_yy_volume_cache_invalidate() bind(C)
  end subroutine
  subroutine cube_yy_wall_cache_metrics(counts) bind(C)
   import c_int64_t
   integer(c_int64_t),intent(out)::counts(4)
  end subroutine
  subroutine cube_yy_wall_invalidate() bind(C)
  end subroutine
  subroutine cube_yy_wall_cache_enable(enabled) bind(C)
   import c_int
   integer(c_int),value::enabled
  end subroutine
  integer(c_int) function cube_yy_exact_unique_points(keys,n,ranks,map,unique) bind(C)
   import c_int,c_int64_t
   integer(c_int64_t),intent(in)::keys(3,*)
   integer(c_int),value::n,ranks
   integer(c_int),intent(out)::map(*)
   integer(c_int64_t),intent(out)::unique(3,*)
  end function
  integer(c_int) function cube_yy_wall_configure(length,omega) bind(C)
   import c_double,c_int
   real(c_double),value::length,omega
  end function
  integer(c_int) function cube_yy_wall_ghost(position,time,source,ghost) bind(C)
   import c_double,c_int
   real(c_double),intent(in)::position(3),source(3)
   real(c_double),value::time
   real(c_double),intent(out)::ghost(3)
  end function
 end interface
contains
 function point_at_held(point)result(p)
  real(c_double),intent(in)::point(3)
  real(c_double)::p(3),a,c,s
  ! Full YY IGRF rotates rigidly about +Z. Native capPhi is held in magnetic
  ! coordinates: rotate the QUERY back to the same held-epoch YY cache.
  a=spin_omega*(query_time-held_time);c=cos(a);s=sin(a)
  p=[c*point(1)+s*point(2),-s*point(1)+c*point(2),point(3)]
 end function
 subroutine yy_boundary_set_time(time_s,status)
  real(c_double),intent(in)::time_s
  integer,intent(out)::status
  status=1
  if(.not.ieee_is_finite(time_s).or.time_s<0) return
  if(installed.and.time_s<held_time) return
  if(time_s/=query_time) then
   query_time=time_s
   call invalidate_point_values()
   call cube_yy_wall_invalidate()
  endif
  status=0
 end subroutine
 subroutine yy_boundary_wall_cache_metrics(counts)
  integer(int64),intent(out)::counts(4)
  call cube_yy_wall_cache_metrics(counts)
 end subroutine
 subroutine yy_boundary_configure(model_length,omega,status)
  real(c_double),intent(in)::model_length,omega
  integer,intent(out)::status
  status=cube_yy_wall_configure(model_length,omega)
  if(status/=0) return
  length_m=model_length;spin_omega=omega;query_time=-1
  call reset_point_cache()
 end subroutine
 subroutine yy_boundary_install(x,y,pn,ps,time_s,status)
  real(c_double),intent(in)::x(:,:),y(:,:),pn(:,:),ps(:,:),time_s
  integer,intent(out)::status
  status=1
  if(length_m<=0.or..not.ieee_is_finite(time_s).or.time_s<0) return
  if(any(shape(x)/=shape(y)).or.any(shape(x)/=shape(pn)).or.any(shape(x)/=shape(ps))) return
  if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(y)).or. &
    .not.all(ieee_is_finite(pn)).or..not.all(ieee_is_finite(ps))) return
  cap_x=x;cap_y=y;north=pn;south=ps;held_time=time_s;query_time=time_s;installed=.true.;status=0
  call exact_cap_configure(cap_x,cap_y,exact_cap_usable)
  call invalidate_point_values()
  call cube_yy_wall_invalidate()
  call cube_yy_volume_cache_invalidate()
 end subroutine
 subroutine yy_boundary_potential_direct(position_m,value,status)
  real(c_double),intent(in)::position_m(3)
  real(c_double),intent(out)::value
  integer,intent(out)::status
  value=0;status=0
  if(.not.installed) return
  status=cube_yy_volume_potential(point_at_held(position_m/length_m),value)
 end subroutine
 subroutine reset_point_cache()
  if(allocated(point_used)) then
   point_used=.false.;point_valid=.false.
  endif
  point_count=0;cache_hits=0;cache_misses=0
 end subroutine
 subroutine invalidate_point_values()
  ! Retain exact query coordinates; every cap installation invalidates values,
  ! even when its time is unchanged. The shape/geometry changes clear all keys.
  if(allocated(point_valid)) point_valid=.false.
  cache_hits=0;cache_misses=0
 end subroutine
 subroutine yy_boundary_collective_refresh(comm,status)
  ! All current MHD owners refresh the exact released volume cache at the
  ! held potential epoch. This derived cache is rebuilt after own restart.
  type(MPI_Comm),intent(in)::comm
  integer,intent(out)::status
  status=cube_yy_volume_cache_refresh(int(comm%MPI_VAL,c_int),held_time)
 end subroutine
 subroutine yy_boundary_cache_enable(enabled)
  logical,intent(in)::enabled
  cache_enabled=enabled;call reset_point_cache()
  call cube_yy_wall_cache_enable(merge(1_c_int,0_c_int,enabled))
 end subroutine
 subroutine yy_boundary_cache_metrics(hits,misses,entries)
  integer(int64),intent(out)::hits,misses
  integer,intent(out)::entries
  hits=cache_hits;misses=cache_misses;entries=point_count
 end subroutine
 integer function point_slot(key) result(slot)
  integer(int64),intent(in)::key(3)
  integer(int64)::mixed
  mixed=ieor(ieor(key(1),ishftc(key(2),21)),ishftc(key(3),43))
  mixed=ieor(mixed,ishft(mixed,-29))
  slot=int(iand(mixed,int(size(point_used)-1,int64)))+1
  do while(point_used(slot))
   if(all(point_key(:,slot)==key)) return
   slot=mod(slot,size(point_used))+1
  enddo
 end function
 subroutine grow_point_cache()
  integer(int64),allocatable::oldkey(:,:)
  real(c_double),allocatable::oldvalue(:)
  logical,allocatable::oldused(:),oldvalid(:)
  integer::n,i,slot
  if(.not.allocated(point_used)) then
   allocate(point_key(3,1024),point_value(1024),point_used(1024),point_valid(1024))
   point_used=.false.;point_valid=.false.;return
  endif
  if(2*point_count<size(point_used)) return
  n=2*size(point_used)
  call move_alloc(point_valid,oldvalid)
  call move_alloc(point_key,oldkey);call move_alloc(point_value,oldvalue);call move_alloc(point_used,oldused)
  allocate(point_key(3,n),point_value(n),point_used(n),point_valid(n))
  point_used=.false.;point_valid=.false.
  do i=1,size(oldused)
   if(.not.oldused(i)) cycle
   slot=point_slot(oldkey(:,i));point_key(:,slot)=oldkey(:,i);point_value(slot)=oldvalue(i);point_used(slot)=.true.
   point_valid(slot)=oldvalid(i)
  enddo
 end subroutine
 subroutine yy_boundary_potential(position_m,value,status)
  real(c_double),intent(in)::position_m(3)
  real(c_double),intent(out)::value
  integer,intent(out)::status
  integer(int64)::key(3)
  integer::slot
  if(.not.cache_enabled.or..not.installed.or..not.all(ieee_is_finite(position_m))) then
   call yy_boundary_potential_direct(position_m,value,status);return
  endif
  call grow_point_cache();key=transfer(position_m,key);slot=point_slot(key)
  if(point_used(slot).and.point_valid(slot)) then
   value=point_value(slot);status=0;cache_hits=cache_hits+1;return
  endif
  call yy_boundary_potential_direct(position_m,value,status);cache_misses=cache_misses+1
  if(status/=0) return
  if(.not.point_used(slot)) point_count=point_count+1
  point_key(:,slot)=key;point_value(slot)=value;point_used(slot)=.true.;point_valid(slot)=.true.
 end subroutine
 integer(c_int) function cube_yy_potential_callback(point,value) bind(C) result(status)
  real(c_double),intent(in)::point(3)
  real(c_double),intent(out)::value
  integer::st
  ! C wall points already use model-length units; avoid an SI round trip.
  value=0;st=0
  if(installed) st=cube_yy_volume_potential(point_at_held(point),value)
  status=st
 end function
 integer(c_int) function cube_yy_footpoint_potential(th,ph,value) bind(C) result(status)
  real(c_double),value::th,ph
  real(c_double),intent(out)::value
  real(c_double)::t,p,w(4),dw(2,4),v(4),fm(3)
  integer::ij(2),st,j,k
  value=0;status=1
  if(.not.installed) return
  ! Native X2ION/Y2ION and N/S tables live in magnetic coordinates.
  ! The released tracer endpoints are geographic; invert the same LFM frame.
  fm=matmul(transpose(yy_native_cap_frame(held_time)), &
   [sin(th)*cos(ph),sin(th)*sin(ph),cos(th)])
  t=atan2(norm2(fm(:2)),abs(fm(3)));p=atan2(fm(2),fm(1))
  if(exact_cap_usable)then
   call exact_cap_weights(t*[cos(p),sin(p)],ij,w,dw,st)
  else
  call mi_cap_weights(cap_x,cap_y,t*[cos(p),sin(p)],ij,w,dw,st)
  endif
  status=st
  if(st==2) then
   status=0;return
  endif
  if(st/=0) return
  j=ij(1);k=ij(2)
  if(fm(3)>0) then
   v=[north(j,k),north(j+1,k),north(j,k+1),north(j+1,k+1)]
  else
   v=[south(j,k),south(j+1,k),south(j,k+1),south(j+1,k+1)]
  endif
  value=dot_product(w,v)
  if(.not.ieee_is_finite(value)) status=1
 end function
 subroutine yy_boundary_volume_domain(low_cm,high_cm,rion_m,status)
  real(c_double),intent(in)::low_cm,high_cm,rion_m
  integer,intent(out)::status
  status=1
  if(length_m<=0) return
  status=cube_yy_volume_cache_configure(max(rion_m/length_m,.98_c_double*.01_c_double*low_cm/length_m), &
   1.02_c_double*.01_c_double*high_cm/length_m,rion_m/length_m)
 end subroutine
 subroutine yy_boundary_ghost(position_cm,time_s,source_cm_s,ghost_cm_s,status)
  real(c_double),intent(in)::position_cm(3),source_cm_s(3),time_s
  real(c_double),intent(out)::ghost_cm_s(3)
  integer,intent(out)::status
  call yy_boundary_set_time(time_s,status);if(status/=0)return
  status=cube_yy_wall_ghost(.01_c_double*position_cm,time_s,.01_c_double*source_cm_s,ghost_cm_s)
  if(status==0) ghost_cm_s=100*ghost_cm_s
 end subroutine
end module
