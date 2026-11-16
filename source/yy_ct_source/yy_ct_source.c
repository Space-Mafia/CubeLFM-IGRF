/* Exact donor CT-flux -> Cartesian-cell-field routine; original Cube grid.
 * This field is an MI source provider; native MHD recovery stays unchanged. */
#include "yy_ct_source.h"
#include "nonorthogonal_step.h"
#include <math.h>
#include <stdlib.h>
void cube_yy_ct_source_destroy(void *handle) {
 cube_yy_ct_source *p=handle;if(!p)return;
 for(int d=0;d<3;++d)free(p->flux[d]);
 free(p->magnetic);gamera_no_grid_destroy(&p->grid);free(p);
}
int cube_yy_ct_source_prepare(int ni,int nj,int nk,double length_m,
                             const double *vertices_cm,void **handle) {
 if(!handle)return 1;
 *handle=NULL;
 if(ni<4||nj<4||nk<4||!vertices_cm||!isfinite(length_m)||length_m<=0.)return 1;
 cube_yy_ct_source *p=calloc(1,sizeof(*p));if(!p)return 2;
 p->n[0]=ni;p->n[1]=nj;p->n[2]=nk;p->length_m=length_m;
 size_t n[3]={(size_t)ni,(size_t)nj,(size_t)nk},v[3]={n[0]+1,n[1]+1,n[2]+1};
 size_t nv=v[0]*v[1]*v[2],nc=n[0]*n[1]*n[2];
 gamera_no_vec3 *xyz=calloc(nv,sizeof(*xyz));p->magnetic=calloc(nc,sizeof(*p->magnetic));
 if(!xyz||!p->magnetic){free(xyz);cube_yy_ct_source_destroy(p);return 2;}
 for(size_t i=0;i<v[0];++i)for(size_t j=0;j<v[1];++j)for(size_t k=0;k<v[2];++k)
  for(int c=0;c<3;++c){double x=vertices_cm[c+3*(i+v[0]*(j+v[1]*k))]*.01;
   if(!isfinite(x)){free(xyz);cube_yy_ct_source_destroy(p);return 3;}
   xyz[gamera_no_index3(v,i,j,k)].value[c]=x;}
 int bad=gamera_no_grid_create(n,xyz,&p->grid);free(xyz);
 if(bad){cube_yy_ct_source_destroy(p);return 4;}
 for(int d=0;d<3;++d){
  size_t count=gamera_no_element_count3(p->grid.face[d].extent);
  p->flux[d]=calloc(count,sizeof(double));
  if(!p->flux[d]){cube_yy_ct_source_destroy(p);return 2;}
 }
 *handle=p;return 0;
}
int cube_yy_ct_source_apply(void *handle,const double *total_phi_cgs,
                            const double *background_phi_cgs,double *cell_gauss) {
 cube_yy_ct_source *p=handle;if(!p||!total_phi_cgs||!background_phi_cgs||!cell_gauss)return 1;
 const size_t ni=p->n[0],nj=p->n[1],nk=p->n[2],vi=ni+1,vj=nj+1,vk=nk+1;
 for(int d=0;d<3;++d){
  const size_t *e=p->grid.face[d].extent;
  for(size_t i=0;i<e[0];++i)for(size_t j=0;j<e[1];++j)for(size_t k=0;k<e[2];++k){
   size_t q[3]={i+1,j+1,k+1};--q[d];
   size_t source=q[0]+vi*(q[1]+vj*(q[2]+vk*d));
   double flux=(total_phi_cgs[source]-background_phi_cgs[source])*1.e-8;
   if(!isfinite(flux))return 2;
   p->flux[d][gamera_no_index3(e,i,j,k)]=flux;
  }
 }
 const double *flux[3]={p->flux[0],p->flux[1],p->flux[2]};
 if(gamera_no_recover_magnetic_field(&p->grid,flux,p->magnetic))return 3;
 for(size_t i=0;i<ni;++i)for(size_t j=0;j<nj;++j)for(size_t k=0;k<nk;++k){
  size_t cell=gamera_no_index3(p->grid.cell_extent,i,j,k);
  for(int c=0;c<3;++c){double b=p->magnetic[cell].value[c]*1.e4;
   if(!isfinite(b))return 4;
   cell_gauss[i+ni*(j+nj*(k+nk*c))]=b;}
 }
 return 0;
}
