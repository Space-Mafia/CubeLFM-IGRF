 subroutine yy_map_global_points(mesh,length_m,global_parallel,t,theta,radii,rion,timefeet,samples,status)
  type(cube_meridian),intent(in)::mesh
  real(real64),intent(in)::length_m,t,theta(:),radii(:),rion,timefeet(:,:)
  real(real64),contiguous,intent(in)::global_parallel(0:,0:,0:,0:,:)
  real(real64),intent(out)::samples(:,:)
  integer,intent(out)::status
  real(real64)::foot(3),point(3),reference(3),logical_xyz(3),direction(3),bmag,imag,z,rad,a,b
  integer::i,st,clip,f,nang,channel
  real(real64),parameter::pi=acos(-1._real64)
  status=1;samples=0;nang=size(global_parallel,2)-2
  do i=1,size(timefeet,2)
   foot=timefeet(:,i)
   call yy_shell_intersection(rion*foot,t,theta,radii,point,st)
   samples(5,i)=real(st,real64)
   if(st/=0)cycle
   call mesh%reference(point*100._real64,rad,reference,st)
   if(st/=0)return
   call cube_inverse(reference,f,a,b,st)
   if(st/=0)return
   logical_xyz=[rad*mesh%ni+.5_real64,(a/(pi/4)+1)*nang/2+.5_real64,(b/(pi/4)+1)*nang/2+.5_real64]
   st=interpolate(logical_xyz,int(f,c_int),int(nang,c_int),global_parallel(:,:,:,:,1),z,clip)
   if(st/=0.or.clip/=0)return
   st=cube_yy_igrf_field(point/length_m,t,direction,bmag)
   if(st/=0)return
   st=cube_yy_igrf_field(rion/length_m*foot,t,direction,imag)
   if(st/=0)return
   samples(3,i)=z*(imag/bmag)
   samples(4,i)=dot_product(direction,foot)
   do channel=1,2
    st=interpolate(logical_xyz,int(f,c_int),int(nang,c_int),global_parallel(:,:,:,:,channel+1),samples(channel,i),clip)
    if(st/=0.or.clip/=0.or.samples(channel,i)<=0)return
   enddo
  enddo
  if(all(ieee_is_finite(samples)))status=0
 end subroutine
