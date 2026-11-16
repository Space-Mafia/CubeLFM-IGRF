/* Original Cube coordinates/data; exact qualified donor geometry and curl.
 * Geometry lives for the model lifetime, NOT one potential solve. */
#include "nonorthogonal_mi.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
typedef struct {
  int nj,nk;
  double length_m;
  gamera_no_grid grid;
  gamera_no_vec3 *magnetic;
} cube_yy_current;
static size_t fv(size_t c,size_t i,size_t j,size_t k,size_t ni,size_t nj) {
 return c+3*(i+ni*(j+nj*k));
}
extern int cube_yy_igrf_field(const double[3],double,double[3],double*);
void cube_yy_current_destroy(void *handle) {
 cube_yy_current *p=handle;if(!p)return;
 gamera_no_grid_destroy(&p->grid);free(p->magnetic);free(p);
}
int cube_yy_current_prepare(int nj,int nk,double length_m,const double *vertices_cm,void **handle) {
 if(!handle)return 1;
 *handle=NULL;
 if(nj<3||nk<3||!vertices_cm||!isfinite(length_m)||length_m<=0.)return 1;
 cube_yy_current *p=calloc(1,sizeof(*p));if(!p)return 2;
 p->nj=nj;p->nk=nk;p->length_m=length_m;
 const size_t n[3]={8,(size_t)nj+4,(size_t)nk+4};
 const size_t v[3]={9,(size_t)nj+5,(size_t)nk+5};
 size_t nv=v[0]*v[1]*v[2],nc=n[0]*n[1]*n[2];
 gamera_no_vec3 *xyz=calloc(nv,sizeof(*xyz));p->magnetic=calloc(nc,sizeof(*p->magnetic));
 if(!xyz||!p->magnetic){free(xyz);cube_yy_current_destroy(p);return 2;}
 for(size_t i=0;i<v[0];++i)for(size_t j=0;j<v[1];++j)for(size_t k=0;k<v[2];++k)
  for(size_t c=0;c<3;++c){double x=.01*vertices_cm[fv(c,i,j,k,v[0],v[1])];
   if(!isfinite(x)){free(xyz);cube_yy_current_destroy(p);return 3;}
   xyz[gamera_no_index3(v,i,j,k)].value[c]=x;}
 int bad=gamera_no_grid_create(n,xyz,&p->grid);free(xyz);
 if(bad){cube_yy_current_destroy(p);return 4;}
 *handle=p;return 0;
}
int cube_yy_current_apply(void *handle,const double *residual_gauss,double time_s,
                          double *current_si,double *parallel_si) {
 cube_yy_current *p=handle;
 if(!p||!residual_gauss||!current_si||!parallel_si||!isfinite(time_s))return 1;
 const size_t n[3]={8,(size_t)p->nj+4,(size_t)p->nk+4};
 for(size_t i=0;i<n[0];++i)for(size_t j=0;j<n[1];++j)for(size_t k=0;k<n[2];++k)
  for(size_t c=0;c<3;++c){double x=1.e-4*residual_gauss[fv(c,i,j,k,n[0],n[1])];
   if(!isfinite(x))return 2;
   p->magnetic[gamera_no_index3(n,i,j,k)].value[c]=x;}
 const double mu0=4.e-7*acos(-1.);
 for(size_t k=1;k<=(size_t)p->nk+2;++k)for(size_t j=1;j<=(size_t)p->nj+2;++j)for(size_t i=1;i<=6;++i){
  gamera_no_vec3 current;
  if(gamera_no_cell_current_from_residual(&p->grid,p->magnetic,i,j,k,&current))return 3;
  size_t cell=gamera_no_index3(p->grid.cell_extent,i,j,k);
  double point[3],direction[3],mag;
  for(int c=0;c<3;++c)point[c]=p->grid.cell[cell].centroid.value[c]/p->length_m;
  if(cube_yy_igrf_field(point,time_s,direction,&mag)||!(mag>0.))return 4;
  double parallel=0.;gamera_no_vec3 true_bg;
  for(int c=0;c<3;++c)true_bg.value[c]=direction[c]*mag;
  for(size_t c=0;c<3;++c){double x=current.value[c]/mu0;
   if(!isfinite(x))return 5;
   current_si[fv(c,i-1,j-1,k-1,6,p->nj+2)]=x;
   parallel+=current.value[c]*true_bg.value[c];
  }
  size_t index=(i-1)+6*((j-1)+(p->nj+2)*(k-1));
  parallel_si[index]=(parallel/mag)/mu0;
  if(!isfinite(parallel_si[index]))return 5;
 }
 return 0;
}
