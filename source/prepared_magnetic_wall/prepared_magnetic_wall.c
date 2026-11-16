/* Equivalent execution of the released magnetic wall. The original object
 * supplies geometry, configuration and the unchanged fallback/oracle. Only
 * external Fortran create/apply/destroy calls are redirected by --wrap.
 * No field interpolation, phase recurrence, changed harmonic sum or mask. */
#include "yy_magnetic_wall.h"
#include <mpi.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GAMERA_MI_COUPLING
#define norm_config cube_yy_ct_norm_config
#define mi_config cube_yy_ct_mi_config
#include "config.h"
#include "background_coeffs.h"
extern Norm_t cube_yy_ct_norm_config;
extern MI_config_t cube_yy_ct_mi_config;
extern int __real_cube_yy_magnetic_wall_create(int,int,int,const double*,const double*,const double*,void**);
extern int __real_cube_yy_magnetic_wall_apply(void*,int,double,double,double,double*);
extern void __real_cube_yy_magnetic_wall_destroy(void*);
typedef struct {
  double phi,ct,st,cp,sp;
  double times[4],flux[4];
  unsigned used,next;
} point_plan;
typedef struct {
  cube_yy_magnetic_wall *original;
  size_t count,stride;
  int degree;
  double signature[5];
  point_plan *points;
  double *spatial,*terms;
  size_t term_stride;int term_count;
  int term_n[SH_N_MAX*(SH_N_MAX+1)/2],term_m[SH_N_MAX*(SH_N_MAX+1)/2];
  double term_g[SH_N_MAX*(SH_N_MAX+1)/2],term_h[SH_N_MAX*(SH_N_MAX+1)/2];
  uint64_t evaluations,hits;
} wall_plan;
_Static_assert(sizeof(MPI_Fint)==sizeof(int),"Existing MPI F08/C handle ABI");
static wall_plan *active_inner_plan,*distributed_plan;
static int distributed_clone,distributed_rank,distributed_size,distributed_count;
static int distributed_counts[4],distributed_offsets[4];
static double *distributed_local,*distributed_current,*distributed_previous;
static double distributed_time,distributed_old_time;
static int distributed_current_ready,distributed_previous_ready,distributed_qualified_epochs;

static void signature(double out[5]) {
  out[0]=norm_config.x_Norm;out[1]=norm_config.B_Norm;
  out[2]=norm_config.Time_Norm;out[3]=mi_config.omega_rad_s;
  out[4]=mi_config.rotation_phase0_rad;
}
static int prepare_point(wall_plan *w,size_t index) {
  const cube_yy_wall_face *f=&w->original->face[index];
  point_plan *p=&w->points[index];
  const int degree=w->degree,width=degree+1;
  double *P=w->spatial+index*w->stride;
  double *dP=P+width*width,*power=dP+width*width;
  const double x=f->point.value[0],y=f->point.value[1],z=f->point.value[2];
  const double r=sqrt(x*x+y*y+z*z);
  if(!(r>=1.e-10)||!isfinite(r))return 1;
  const double theta=acos(z/r);
  p->phi=atan2(y,x);p->ct=cos(theta);p->st=sin(theta);
  if(p->st<1.e-10)p->st=1.e-10;
  p->cp=cos(p->phi);p->sp=sin(p->phi);
  P[0]=1.;dP[0]=0.;
  /* Original donor recurrence, operations and n/m ordering. */
  for(int n=1;n<=degree;++n)for(int m=0;m<=n;++m){
    const size_t q=n*width+m,old=(n-1)*width+m;
    if(n==m){
      const double S=(n==1)?1.:sqrt((double)(2*n-1)/(2.*n));
      P[q]=S*p->st*P[(n-1)*width+n-1];
      dP[q]=S*(p->ct*P[(n-1)*width+n-1]+p->st*dP[(n-1)*width+n-1]);
    }else{
      const double a=(2.*n-1.)/sqrt((double)(n*n-m*m));
      const double b=sqrt((double)((n-1)*(n-1)-m*m)/(double)(n*n-m*m));
      P[q]=a*p->ct*P[old]-b*(n>1?P[(n-2)*width+m]:0.);
      dP[q]=a*(-p->st*P[old]+p->ct*dP[old])-b*(n>1?dP[(n-2)*width+m]:0.);
    }
  }
  const double R_sim=(R_PLANET_KM*1000.)/norm_config.x_Norm;
  for(int n=1;n<=degree;++n){const double ratio=R_sim/r;power[n]=pow(ratio,n+2);}
  double *t=w->terms+index*w->term_stride;
  for(int z=0;z<w->term_count;++z){int n=w->term_n[z],m=w->term_m[z];
    t[3*z]=((n+1)*power[n])*P[n*width+m];
    t[3*z+1]=power[n]*dP[n*width+m];
    t[3*z+2]=(power[n]*P[n*width+m])*m;
  }
  return 0;
}
static int prepared_flux(wall_plan *w,size_t index,double time_s,double *flux) {
  point_plan *p=&w->points[index];
  for(unsigned i=0;i<p->used;++i)if(p->times[i]==time_s){
    *flux=p->flux[i];++w->hits;return 0;
  }
  const int degree=w->degree;
  const double time=time_s/w->signature[2];
  const double omega_code=w->signature[3]*w->signature[2];
  const double phi_eff=p->phi-(w->signature[4]+omega_code*time);
  const double nT_to_code=1.e-9/w->signature[1];
  double Br_n=0.,Bth_n=0.,Bph_n=0.;
  double cos_order[SH_N_MAX+1],sin_order[SH_N_MAX+1];
  unsigned char order_ready[SH_N_MAX+1]={0};
  const double *terms=w->terms+index*w->term_stride;
  for(int z=0;z<w->term_count;++z){
    const int m=w->term_m[z];const double g=w->term_g[z],h=w->term_h[z];
    if(!order_ready[m]){cos_order[m]=cos(m*phi_eff);sin_order[m]=sin(m*phi_eff);order_ready[m]=1;}
    const double cos_mp=cos_order[m],sin_mp=sin_order[m];
    Br_n+=terms[3*z]*(g*cos_mp+h*sin_mp);
    Bth_n-=terms[3*z+1]*(g*cos_mp+h*sin_mp);
    Bph_n+=terms[3*z+2]*(g*sin_mp-h*cos_mp)/p->st;
  }
  double b[3];
  b[0]=(Br_n*p->st*p->cp+Bth_n*p->ct*p->cp-Bph_n*p->sp)*nT_to_code;
  b[1]=(Br_n*p->st*p->sp+Bth_n*p->ct*p->sp+Bph_n*p->cp)*nT_to_code;
  b[2]=(Br_n*p->ct-Bth_n*p->st)*nT_to_code;
  for(int a=0;a<3;++a){b[a]=b[a]*w->signature[1]*1.e4;if(!isfinite(b[a]))return 1;}
  const cube_yy_wall_face *f=&w->original->face[index];
  *flux=b[0]*f->vector[0]+b[1]*f->vector[1]+b[2]*f->vector[2];
  if(!isfinite(*flux))return 1;
  const unsigned slot=p->next;p->times[slot]=time_s;p->flux[slot]=*flux;
  p->next=(slot+1)&3u;if(p->used<4)++p->used;++w->evaluations;
  return 0;
}
static void touch_local_pages(void *memory,size_t bytes) {
 if(!bytes)return;volatile unsigned char *p=memory;
 for(size_t i=0;i<bytes;i+=4096)p[i]=0;p[bytes-1]=0;
}
static int worker_status(MPI_Comm comm,int local) {
 int global=0;return MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MAX,comm)==MPI_SUCCESS?global:91;
}
int cube_yy_distributed_wall_prepare(int cf) {
 MPI_Comm comm=MPI_Comm_f2c((MPI_Fint)cf);if(comm==MPI_COMM_NULL||distributed_plan)return 1;
 if(MPI_Comm_rank(comm,&distributed_rank)!=MPI_SUCCESS||MPI_Comm_size(comm,&distributed_size)!=MPI_SUCCESS)return 91;
 if(distributed_size!=4)return 1;
 int header[8]={0},local=0;double values[6]={0};
 if(distributed_rank==0){
  wall_plan *p=active_inner_plan;
  if(!p||!p->points||!p->terms||p->count>INT_MAX/sizeof(point_plan)||p->count>INT_MAX/sizeof(cube_yy_wall_face)||p->count*p->term_stride>INT_MAX)local=1;
  else{
   header[0]=p->degree;header[1]=p->term_count;header[2]=p->original->ng;header[3]=p->original->nj;header[4]=p->original->nk;header[5]=(int)p->count;header[6]=(int)sizeof(point_plan);header[7]=(int)sizeof(cube_yy_wall_face);
   memcpy(values,p->signature,5*sizeof(double));values[5]=p->original->length_cm;
  }
 }
 if(MPI_Bcast(&local,1,MPI_INT,0,comm)!=MPI_SUCCESS)return 91;if(local)return local;
 if(MPI_Bcast(header,8,MPI_INT,0,comm)!=MPI_SUCCESS||MPI_Bcast(values,6,MPI_DOUBLE,0,comm)!=MPI_SUCCESS)return 91;
 if(header[6]!=(int)sizeof(point_plan)||header[7]!=(int)sizeof(cube_yy_wall_face)||header[1]<0||header[1]>SH_N_MAX*(SH_N_MAX+1)/2||header[5]<1)local=1;
 if(worker_status(comm,local))return 1;
 distributed_clone=distributed_rank!=0;distributed_plan=distributed_clone?calloc(1,sizeof(wall_plan)):active_inner_plan;
 if(!distributed_plan)local=1;
 if(!local&&distributed_clone){
  wall_plan *p=distributed_plan;p->degree=header[0];p->term_count=header[1];p->count=(size_t)header[5];p->term_stride=(size_t)3*p->term_count;memcpy(p->signature,values,5*sizeof(double));
  p->original=calloc(1,sizeof(cube_yy_magnetic_wall));
  p->points=malloc(p->count*sizeof(point_plan));p->terms=malloc((p->count*p->term_stride?p->count*p->term_stride:1)*sizeof(double));
  if(!p->original||!p->points||!p->terms)local=1;
  if(p->original){p->original->ng=header[2];p->original->nj=header[3];p->original->nk=header[4];p->original->length_cm=values[5];p->original->face=malloc(p->count*sizeof(cube_yy_wall_face));if(!p->original->face)local=1;}
 }
 distributed_count=header[5];
 for(int i=0;i<4;++i){distributed_offsets[i]=distributed_count*i/4;distributed_counts[i]=distributed_count*(i+1)/4-distributed_offsets[i];}
 distributed_local=malloc((size_t)distributed_counts[distributed_rank]*sizeof(double));if(!distributed_local)local=1;
 if(distributed_rank==0){distributed_current=malloc((size_t)distributed_count*sizeof(double));distributed_previous=malloc((size_t)distributed_count*sizeof(double));if(!distributed_current||!distributed_previous)local=1;}
 if(worker_status(comm,local))return 1;
 wall_plan *p=distributed_plan;
 if(distributed_clone){
  touch_local_pages(p->points,p->count*sizeof(point_plan));touch_local_pages(p->terms,p->count*p->term_stride*sizeof(double));touch_local_pages(p->original->face,p->count*sizeof(cube_yy_wall_face));
 }
 if(MPI_Bcast(p->term_n,p->term_count,MPI_INT,0,comm)!=MPI_SUCCESS||MPI_Bcast(p->term_m,p->term_count,MPI_INT,0,comm)!=MPI_SUCCESS||MPI_Bcast(p->term_g,p->term_count,MPI_DOUBLE,0,comm)!=MPI_SUCCESS||MPI_Bcast(p->term_h,p->term_count,MPI_DOUBLE,0,comm)!=MPI_SUCCESS)return 91;
 if(MPI_Bcast(p->points,(int)(p->count*sizeof(point_plan)),MPI_BYTE,0,comm)!=MPI_SUCCESS||MPI_Bcast(p->terms,(int)(p->count*p->term_stride),MPI_DOUBLE,0,comm)!=MPI_SUCCESS||MPI_Bcast(p->original->face,(int)(p->count*sizeof(cube_yy_wall_face)),MPI_BYTE,0,comm)!=MPI_SUCCESS)return 91;
 distributed_current_ready=distributed_previous_ready=distributed_qualified_epochs=0;
 return 0;
}
static int refresh_distributed_epoch(MPI_Comm comm,double time,double *root_output) {
 int local=0,first=distributed_offsets[distributed_rank],n=distributed_counts[distributed_rank];
 for(int i=0;i<n;++i)if(prepared_flux(distributed_plan,(size_t)(first+i),time,distributed_local+i)){local=1;break;}
 int status=worker_status(comm,local);if(status)return status;
 return MPI_Gatherv(distributed_local,n,MPI_DOUBLE,root_output,distributed_counts,distributed_offsets,MPI_DOUBLE,0,comm)==MPI_SUCCESS?0:91;
}
int cube_yy_distributed_wall_prefetch(int cf,double time,double dt_old) {
 MPI_Comm comm=MPI_Comm_f2c((MPI_Fint)cf);if(comm==MPI_COMM_NULL||!distributed_plan)return 1;
 double clocks[2]={time,dt_old},root_clocks[2]={time,dt_old};
 if(MPI_Bcast(root_clocks,2,MPI_DOUBLE,0,comm)!=MPI_SUCCESS)return 91;
 int local=memcmp(clocks,root_clocks,sizeof(clocks))||!isfinite(time)||!isfinite(dt_old),status=worker_status(comm,local);if(status)return status;
 distributed_current_ready=distributed_previous_ready=0;
 status=refresh_distributed_epoch(comm,time,distributed_current);if(status)return status;
 distributed_time=time;distributed_current_ready=1;
 if(dt_old>0.){
  status=refresh_distributed_epoch(comm,time-dt_old,distributed_previous);if(status)return status;
  distributed_old_time=time-dt_old;distributed_previous_ready=1;
 }
 return 0;
}
int cube_yy_distributed_wall_check(int cf,double time) {
 MPI_Comm comm=MPI_Comm_f2c((MPI_Fint)cf);if(comm==MPI_COMM_NULL||!distributed_plan)return 1;
 int local=0;
 if(distributed_rank==0){
  if(!distributed_current_ready||distributed_time!=time)local=1;
  for(int i=0;i<distributed_count&&!local;++i){double value;
   if(prepared_flux(distributed_plan,(size_t)i,time,&value)||memcmp(&value,distributed_current+i,sizeof(double)))local=2;
  }
 }
 if(MPI_Bcast(&local,1,MPI_INT,0,comm)!=MPI_SUCCESS)return 91;if(local)return local;
 ++distributed_qualified_epochs;
 if(distributed_qualified_epochs==8&&distributed_rank==0){
  int worldrank=-1;if(MPI_Comm_rank(MPI_COMM_WORLD,&worldrank)!=MPI_SUCCESS)local=91;char name[128];snprintf(name,sizeof(name),"original-wall-worker-word-r%03d.json",worldrank);
  FILE *file=local?NULL:fopen(name,"wx");if(!file&&local==0)local=3;
  if(file){
   if(fprintf(file,"{\"status\":\"ORIGINAL_WALL_ALL_WORKER_WORD_PASS\",\"rank\":%d,\"epochs\":8,\"values\":%d,\"immutable_geometry_or_tracer_recomputed\":false}\n",worldrank,distributed_count*8)<0)local=3;
   if(fclose(file))local=3;
  }
 }
 if(MPI_Bcast(&local,1,MPI_INT,0,comm)!=MPI_SUCCESS)return 91;return local;
}
void cube_yy_distributed_wall_release(void) {
 if(distributed_clone&&distributed_plan){
  free(distributed_plan->points);free(distributed_plan->terms);
  if(distributed_plan->original){free(distributed_plan->original->face);free(distributed_plan->original);}
  free(distributed_plan);
 }
 distributed_plan=NULL;distributed_clone=0;
 free(distributed_local);free(distributed_current);free(distributed_previous);
 distributed_local=distributed_current=distributed_previous=NULL;distributed_current_ready=distributed_previous_ready=distributed_qualified_epochs=0;
}
void __wrap_cube_yy_magnetic_wall_destroy(void *handle) {
  wall_plan *p=handle;if(!p)return;
  if(active_inner_plan==p)active_inner_plan=NULL;
  fprintf(stdout,"FOUR_SLOT_EXACT_WALL_COUNTS points=%zu evaluations=%llu hits=%llu\n",p->count,(unsigned long long)p->evaluations,(unsigned long long)p->hits);
  __real_cube_yy_magnetic_wall_destroy(p->original);
  free(p->points);free(p->spatial);free(p->terms);free(p);
}
int __wrap_cube_yy_magnetic_wall_create(int nj,int nk,int ng,const double *vertices,
    const double *vectors,const double *norms,void **handle) {
  if(!handle)return 1;*handle=NULL;
  wall_plan *p=calloc(1,sizeof(*p));if(!p)return 2;
  int st=__real_cube_yy_magnetic_wall_create(nj,nk,ng,vertices,vectors,norms,(void**)&p->original);
  if(st){free(p);return st;}
  p->degree=sh_active_degree();signature(p->signature);
  p->count=(size_t)(ng+1)*(nj+2*ng)*(nk+2*ng)*3;
  p->stride=(size_t)2*(p->degree+1)*(p->degree+1)+p->degree+1;
  p->points=calloc(p->count,sizeof(*p->points));p->spatial=calloc(p->count*p->stride,sizeof(double));
  for(int n=1;n<=p->degree;++n)for(int m=1;m<=n;++m){
    double g=sh_coeffs_g[n][m],h=sh_coeffs_h[n][m];if(g==0.&&h==0.)continue;
    int z=p->term_count++;p->term_n[z]=n;p->term_m[z]=m;p->term_g[z]=g;p->term_h[z]=h;
  }
  p->term_stride=(size_t)3*p->term_count;
  p->terms=calloc(p->count*(p->term_stride?p->term_stride:1),sizeof(double));

  /* Allocation/unsupported degeneracy uses the original operator unchanged. */
  if(!p->points||!p->spatial||!p->terms) {free(p->points);free(p->spatial);free(p->terms);p->points=NULL;p->spatial=NULL;p->terms=NULL;}
  if(p->points)for(size_t i=0;i<p->count;++i)if(prepare_point(p,i)){
    free(p->points);free(p->spatial);free(p->terms);p->points=NULL;p->spatial=NULL;p->terms=NULL;break;
  }
  if(p->points){free(p->spatial);p->spatial=NULL;}
  active_inner_plan=p;*handle=p;return 0;
}
int __wrap_cube_yy_magnetic_wall_apply(void *handle,int ni,double time,double dt_old,
    double predictor_factor,double *bf) {
  wall_plan *p=handle;
  if(!p)return 1;
  const cube_yy_magnetic_wall *w=p->original;
  double current_signature[5];signature(current_signature);
  if(!p->points||memcmp(current_signature,p->signature,sizeof(current_signature))||p->degree!=sh_active_degree())
    return __real_cube_yy_magnetic_wall_apply(w? p->original:NULL,ni,time,dt_old,predictor_factor,bf);
  if(!w||!bf||ni<1||!isfinite(time)||!isfinite(dt_old)||
     !isfinite(predictor_factor)||predictor_factor<0.||
     (predictor_factor>0.&&dt_old<=0.)||100.*norm_config.x_Norm!=w->length_cm)return 1;
  const int ni_full=ni+2*w->ng,nj_full=w->nj+2*w->ng,nk_full=w->nk+2*w->ng;
  for(int d=0;d<3;++d)for(int k=0;k<nk_full;++k)for(int j=0;j<nj_full;++j){
    const size_t si=cube_yy_wall_index(w,w->ng,j,k,d);
    const cube_yy_wall_face *src=&w->face[si];
    const size_t source=w->ng+(size_t)ni_full*(j+nj_full*(k+nk_full*d));
    double src_planet,src_old=0.;
    if(p==distributed_plan&&distributed_current_ready&&distributed_time==time)src_planet=distributed_current[si];
    else if(prepared_flux(p,si,time,&src_planet))return 2;
    if(predictor_factor>0.){
      if(p==distributed_plan&&distributed_previous_ready&&distributed_old_time==time-dt_old)src_old=distributed_previous[si];
      else if(prepared_flux(p,si,time-dt_old,&src_old))return 2;
      src_planet=(1.+predictor_factor)*src_planet-predictor_factor*src_old;
    }
    const double perturbation=bf[source]*src->norm-src_planet;
    if(!isfinite(perturbation))return 2;
    for(int i=0;i<w->ng;++i){
      const size_t gi=cube_yy_wall_index(w,i,j,k,d);
      const cube_yy_wall_face *ghost=&w->face[gi];
      double ghost_planet,ghost_old=0.;
      if(p==distributed_plan&&distributed_current_ready&&distributed_time==time)ghost_planet=distributed_current[gi];
      else if(prepared_flux(p,gi,time,&ghost_planet))return 2;
      if(predictor_factor>0.){
        if(p==distributed_plan&&distributed_previous_ready&&distributed_old_time==time-dt_old)ghost_old=distributed_previous[gi];
        else if(prepared_flux(p,gi,time-dt_old,&ghost_old))return 2;
        ghost_planet=(1.+predictor_factor)*ghost_planet-predictor_factor*ghost_old;
      }
      const double ratio=ghost->scalar/src->scalar;
      const size_t target=i+(size_t)ni_full*(j+nj_full*(k+nk_full*d));
      bf[target]=(perturbation*ratio+ghost_planet)/ghost->norm;
      if(!isfinite(bf[target]))return 2;
    }
  }
  return 0;
}
/* Component metadata only; never affects operator arithmetic. */
int cube_yy_prepared_wall_metrics(void *handle,uint64_t out[4]) {
  wall_plan *p=handle;if(!p||!out)return 1;
  out[0]=p->evaluations;out[1]=p->hits;out[2]=p->count;out[3]=p->points!=NULL;return 0;
}
cube_yy_magnetic_wall *cube_yy_prepared_wall_original(void *handle) {
  wall_plan *p=handle;return p?p->original:NULL;
}
