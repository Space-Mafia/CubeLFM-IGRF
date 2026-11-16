/* Storage permutation and original subtraction/projection expressions.
 * Source fields are SoA; output plane's first axis is contiguous. All
 * original finite/positive guards execute before these private routines.
 * Only the qualified local16x12x12/g4 layout enters this fast path. */
#include <stddef.h>
static int metadata(const int n[3],int g,int d,int kk,int fixed_axis){return n&&n[0]==16&&n[1]==12&&n[2]==12&&g==4&&d>=1&&d<=3&&fixed_axis>=0&&fixed_axis<3&&kk>=1&&kk<=n[fixed_axis];}
int cube_yy_pack_hydro_plane(const int n[3],int g,int direction,int kk,
 const double *restrict q,const double *restrict p,const double *restrict shift,
 double *restrict qb,double *restrict qp){
 if(!metadata(n,g,direction,kk,(direction+1)%3)||!q||!p||!shift||!qb||!qp)return 1;
 int d=direction-1,a=(d+1)%3,b=(d+2)%3,nx[3]={n[0]+2*g,n[1]+2*g,n[2]+2*g};
 size_t stride[3]={1,(size_t)nx[0],(size_t)nx[0]*nx[1]},volume=(size_t)nx[0]*nx[1]*nx[2],plane=(size_t)nx[d]*nx[a];
 size_t start=(size_t)(kk+g-1)*stride[b];
 for(int z=0;z<8;++z)for(int j=0;j<nx[a];++j){
  size_t source=start+(size_t)j*stride[a],out=(size_t)j*nx[d]+(size_t)z*plane;
  for(int i=0;i<nx[d];++i){size_t at=source+(size_t)i*stride[d];
   if(z<5){qb[out+i]=q[at+(size_t)z*volume];qp[out+i]=p[at+(size_t)z*volume];}
   else{double ref=shift[at*3+(size_t)z-5];qb[out+i]=q[at+(size_t)z*volume]-ref;qp[out+i]=p[at+(size_t)z*volume]-ref;}
  }
 }
 return 0;
}
int cube_yy_pack_bfield_plane(const int n[3],int g,int direction,int kk,
 const double *restrict q,const double *restrict p,const double *restrict bf,const double *restrict bp,
 const double *restrict face_ref,const double *restrict area,const double *restrict magnitude,
 const double *restrict uniform,const double *restrict emf,
 double *restrict qb,double *restrict qp,double *restrict bb,double *restrict pp){
 if(!metadata(n,g,direction,kk,direction-1)||!q||!p||!bf||!bp||!face_ref||!area||!magnitude||!uniform||!emf||!qb||!qp||!bb||!pp)return 1;
 int d=direction-1,a=(d+1)%3,b=(d+2)%3,nx[3]={n[0]+2*g,n[1]+2*g,n[2]+2*g};
 size_t stride[3]={1,(size_t)nx[0],(size_t)nx[0]*nx[1]},volume=(size_t)nx[0]*nx[1]*nx[2],plane=(size_t)nx[a]*nx[b];
 size_t start=(size_t)(kk+g-1)*stride[d];
 for(int z=0;z<8;++z)for(int j=0;j<nx[b];++j){
  size_t source=start+(size_t)j*stride[b],out=(size_t)j*nx[a]+(size_t)z*plane;
  for(int i=0;i<nx[a];++i){size_t at=source+(size_t)i*stride[a];qb[out+i]=q[at+(size_t)z*volume];qp[out+i]=p[at+(size_t)z*volume];}
 }
 double r[3];for(int z=0;z<3;++z)r[z]=uniform[z]+emf[z];
 for(int z=0;z<2;++z){int axis=z==0?a:b;
  for(int j=0;j<nx[b];++j){size_t source=start+(size_t)j*stride[b],out=(size_t)j*nx[a]+(size_t)z*plane;
   for(int i=0;i<nx[a];++i){size_t at=source+(size_t)i*stride[a],offset=at+(size_t)axis*volume;
    double dot=0.;for(int l=0;l<3;++l)dot=dot+r[l]*area[3*offset+l];
    double ref=face_ref[offset]+dot/magnitude[offset];bb[out+i]=bf[offset]-ref;pp[out+i]=bp[offset]-ref;
   }
  }
 }
 return 0;
}
