/* Layout/units adapter for the RELEASED rotating magnetic inner wall.
 * No active CT face, mesh, fluid boundary, MI solver or clock is changed.
 * The tagged field bodies and geometry library are reused unchanged. */
#include "yy_magnetic_wall.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define GAMERA_MI_COUPLING
#define GAMERA_PERF_HARMONIC_TRIG 1
#define norm_config cube_yy_ct_norm_config
#define mi_config cube_yy_ct_mi_config
#include "config.h"
#include "background_coeffs.h"
#include "harmonic_geometry_cache.h"
extern Norm_t cube_yy_ct_norm_config;
extern MI_config_t cube_yy_ct_mi_config;
#include "tagged_rotating_harmonics.inc"

int cube_yy_magnetic_wall_field(const double point[3],double time_s,double b_gauss[3]) {
  if(!point||!b_gauss||!isfinite(time_s)||!(norm_config.x_Norm>0.)||
     !(norm_config.B_Norm>0.)||!(norm_config.Time_Norm>0.)) return 1;
  gamera_no_vec3 p,b;
  for(int a=0;a<3;++a){if(!isfinite(point[a]))return 1;p.value[a]=point[a];}
  if(!(p.value[0]*p.value[0]+p.value[1]*p.value[1]+p.value[2]*p.value[2]>0.))return 1;
  evaluate_rotating_harmonics(p,time_s/norm_config.Time_Norm,&b);
  for(int a=0;a<3;++a){b_gauss[a]=b.value[a]*norm_config.B_Norm*1.e4;if(!isfinite(b_gauss[a]))return 1;}
  return 0;
}

void cube_yy_magnetic_wall_destroy(void *handle) {
  cube_yy_magnetic_wall *w=handle;
  if(w){free(w->face);free(w);}
}

int cube_yy_magnetic_wall_create(int nj,int nk,int ng,const double *vertices_cm,
    const double *vectors_cm2,const double *norms_cm2,void **handle) {
  if(!handle)return 1;
  *handle=NULL;
  if(nj<1||nk<1||ng<1||!vertices_cm||!vectors_cm2||!norms_cm2||
      !(norm_config.x_Norm>0.)||!isfinite(mi_config.omega_rad_s)||mi_config.omega_rad_s<=0.)return 1;
  cube_yy_magnetic_wall *w=calloc(1,sizeof(*w));if(!w)return 2;
  w->ng=ng;w->nj=nj;w->nk=nk;w->length_cm=100.*norm_config.x_Norm;
  size_t count=(size_t)(ng+1)*(nj+2*ng)*(nk+2*ng)*3;
  w->face=calloc(count,sizeof(*w->face));if(!w->face){cube_yy_magnetic_wall_destroy(w);return 2;}
  const int nv_i=ng+2,nv_j=nj+2*ng+1;
  /* The donor corner order differs from Fortran bit order. */
  const int bit[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                     {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
  for(int k=0;k<nk+2*ng;++k)for(int j=0;j<nj+2*ng;++j)for(int i=0;i<=ng;++i){
    gamera_no_vec3 v[8];gamera_no_cell_geometry geometry;
    for(int h=0;h<8;++h)for(int a=0;a<3;++a){
      size_t p=a+3*((i+bit[h][0])+nv_i*((j+bit[h][1])+nv_j*(k+bit[h][2])));
      v[h].value[a]=vertices_cm[p]/w->length_cm;
    }
    if(gamera_no_compute_cell_geometry(v,&geometry)){cube_yy_magnetic_wall_destroy(w);return 3;}
    for(int d=0;d<3;++d){
      size_t p=cube_yy_wall_index(w,i,j,k,d);cube_yy_wall_face *f=&w->face[p];
      f->point=geometry.face[d][GAMERA_NO_LOWER].centroid;
      f->scalar=geometry.face[d][GAMERA_NO_LOWER].area*w->length_cm*w->length_cm;
      f->norm=norms_cm2[p];
      for(int a=0;a<3;++a){f->vector[a]=vectors_cm2[a+3*p];
        if(!isfinite(f->vector[a])){cube_yy_magnetic_wall_destroy(w);return 4;}}
      if(!isfinite(f->scalar)||!isfinite(f->norm)||f->scalar<=0.||f->norm<=0.){
        cube_yy_magnetic_wall_destroy(w);return 4;}
    }
  }
  *handle=w;return 0;
}

static int planet_flux(const cube_yy_wall_face *f,double time,double *flux) {
  double b[3];if(cube_yy_magnetic_wall_field(f->point.value,time,b))return 1;
  *flux=b[0]*f->vector[0]+b[1]*f->vector[1]+b[2]*f->vector[2];
  return !isfinite(*flux);
}

int cube_yy_magnetic_wall_apply(void *handle,int ni,double time,double dt_old,
                               double predictor_factor,double *bf) {
  const cube_yy_magnetic_wall *w=handle;
  if(!w||!bf||ni<1||!isfinite(time)||!isfinite(dt_old)||
      !isfinite(predictor_factor)||predictor_factor<0.||
      (predictor_factor>0.&&dt_old<=0.)||100.*norm_config.x_Norm!=w->length_cm)return 1;
  const int ni_full=ni+2*w->ng,nj_full=w->nj+2*w->ng,nk_full=w->nk+2*w->ng;
  for(int d=0;d<3;++d)for(int k=0;k<nk_full;++k)for(int j=0;j<nj_full;++j){
    const cube_yy_wall_face *src=&w->face[cube_yy_wall_index(w,w->ng,j,k,d)];
    size_t source=w->ng+(size_t)ni_full*(j+nj_full*(k+nk_full*d));
    double src_planet,src_old=0.;if(planet_flux(src,time,&src_planet))return 2;
    if(predictor_factor>0.){
      if(planet_flux(src,time-dt_old,&src_old))return 2;
      src_planet=(1.+predictor_factor)*src_planet-predictor_factor*src_old;
    }
    double perturbation=bf[source]*src->norm-src_planet;
    if(!isfinite(perturbation))return 2;
    for(int i=0;i<w->ng;++i){
      const cube_yy_wall_face *ghost=&w->face[cube_yy_wall_index(w,i,j,k,d)];
      double ghost_planet,ghost_old=0.;if(planet_flux(ghost,time,&ghost_planet))return 2;
      if(predictor_factor>0.){
        if(planet_flux(ghost,time-dt_old,&ghost_old))return 2;
        ghost_planet=(1.+predictor_factor)*ghost_planet-predictor_factor*ghost_old;
      }
      double ratio=ghost->scalar/src->scalar;
      size_t target=i+(size_t)ni_full*(j+nj_full*(k+nk_full*d));
      bf[target]=(perturbation*ratio+ghost_planet)/ghost->norm;
      if(!isfinite(bf[target]))return 2;
    }
  }
  return 0;
}
