/* Storage adapter only. The existing released YY routine owns all decisions. */
#include "nonorthogonal_step.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
static double *packed;
static size_t capacity, calls, faces;
static double massmax, energymax;
static size_t offset(size_t ni,size_t nj,size_t nk,size_t j,size_t k,size_t c) {
 return (ni+1U)*(j+(nj+1U)*(k+(nk+1U)*3U*c));
}
int cube_yy_inner_fluid_flux_boundary(int ni,int nj,int nk,double *f) {
 if(ni<1||nj<1||nk<1||f==NULL) return -1;
 size_t nr=(size_t)ni,nt=(size_t)nj,np=(size_t)nk;
 if(nt>SIZE_MAX/np/5U/sizeof(double)) return -1;
 size_t needed=nt*np*5U;
 if(needed>capacity) {
  double *p=realloc(packed,needed*sizeof(double));
  if(!p) return -1;
  packed=p;capacity=needed;
 }
 for(size_t j=1;j<=nt;++j) for(size_t k=1;k<=np;++k) {
  size_t cell=((j-1)*np+k-1)*5U;
  for(size_t c=0;c<5;++c) {
   double v=f[offset(nr,nt,np,j,k,c)];
   if(!isfinite(v)) return -1;
   packed[cell+c]=v;
  }
  if(packed[cell]>0.0) {
   massmax=fmax(massmax,packed[cell]);
   energymax=fmax(energymax,packed[cell+4]);
  }
 }
 gamera_no_storage storage;memset(&storage,0,sizeof(storage));
 storage.face_extent[0][0]=1;storage.face_extent[0][1]=nt;storage.face_extent[0][2]=np;
 storage.fluid_face_flux[0]=packed;
 const size_t lower[3]={0,0,0},upper[3]={1,nt,np};size_t count=0;
 int status=gamera_no_trap_inner_outward_mass_energy_flux(&storage,lower,upper,&count);
 if(status!=0) return status;
 for(size_t j=1;j<=nt;++j) for(size_t k=1;k<=np;++k) {
  size_t cell=((j-1)*np+k-1)*5U;
  f[offset(nr,nt,np,j,k,0)]=packed[cell];
  f[offset(nr,nt,np,j,k,4)]=packed[cell+4];
 }
 ++calls;faces+=count;return 0;
}
void cube_yy_inner_flux_statistics(size_t *ncall,size_t *nface,double *mmax,double *emax) {
 *ncall=calls;*nface=faces;*mmax=massmax;*emax=energymax;
}
