/* Production lifecycle adapter around the released YY volume return operator.
 * Interpolation is byte exact; refresh only substitutes the current Cube MPI
 * communicator for MPI_COMM_WORLD. The field/tracer and C gas wall are reused. */
#include <mpi.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#define CACHE_NR 8
#define CACHE_NTH 91
#define CACHE_NPH 180
#define PI (acos(-1.0))
#define log_info(...) do { fprintf(stdout,__VA_ARGS__); fputc('\n',stdout); } while(0)
#define GAMERA_NO_DIM 3
typedef struct { double value[3]; } gamera_no_vec3;
static int rank;
static MPI_Comm cache_comm=MPI_COMM_NULL;
static struct { double ionosphere_radius_re; } mi_config;
static struct { int solution_ready; double *pot_cache_th_foot,*pot_cache_ph_foot,*pot_cache_value;
 double pot_cache_r_min,pot_cache_r_max; } mi_state;
static double vector_norm(gamera_no_vec3 v) {
 double sum=0.;for(int i=0;i<3;++i)sum+=v.value[i]*v.value[i];return sqrt(fmax(0.,sum));
}
_Static_assert(sizeof(MPI_Fint)==sizeof(int),"Cube F08 communicator ABI requires MPI_Fint/int equality");
static int interpolate_potential(gamera_no_vec3 point, double *potential_v) {
  if (potential_v == NULL || !mi_state.solution_ready) return -1;
  
  // 1. Calculate angles using the TRUE unclamped radius
  double true_radius = vector_norm(point);
  if (!(true_radius > 0.0)) {
      *potential_v = 0.0;
      return 0;
  }
  
  // Bound the ratio to [-1, 1] before passing to acos
  double th = acos(fmax(-1.0, fmin(1.0, point.value[2] / true_radius)));
  double ph = atan2(point.value[1], point.value[0]);
  if (ph < 0.0) ph += 2.0 * PI;

  // 2. Clamp the radius for the cache lookup
  double radius = true_radius;
  if (radius < mi_state.pot_cache_r_min) radius = mi_state.pot_cache_r_min;
  if (radius > mi_state.pot_cache_r_max) radius = mi_state.pot_cache_r_max;
  
  // ... (proceed with interpolation logic) ...

  double dr = (mi_state.pot_cache_r_max - mi_state.pot_cache_r_min) / (double)(CACHE_NR - 1);
  double dth = PI / (double)(CACHE_NTH - 1);
  double dph = (2.0 * PI) / (double)CACHE_NPH;

  double r_log = (radius - mi_state.pot_cache_r_min) / dr;
  double th_log = th / dth;
  double ph_log = ph / dph;

  int i0 = (int)floor(r_log);  if (i0 >= CACHE_NR - 1) i0 = CACHE_NR - 2;
  int j0 = (int)floor(th_log); if (j0 >= CACHE_NTH - 1) j0 = CACHE_NTH - 2;
  int k0 = (int)floor(ph_log) % CACHE_NPH;
  int k1 = (k0 + 1) % CACHE_NPH;

  double rd = r_log - i0;
  double td = th_log - j0;
  double pd = ph_log - floor(ph_log);

  #define GET_CACHE(ri, ti, pi) mi_state.pot_cache_value[((ri) * CACHE_NTH * CACHE_NPH) + ((ti) * CACHE_NPH) + (pi)]

  double c000 = GET_CACHE(i0, j0, k0);   double c100 = GET_CACHE(i0+1, j0, k0);
  double c001 = GET_CACHE(i0, j0, k1);   double c101 = GET_CACHE(i0+1, j0, k1);
  double c010 = GET_CACHE(i0, j0+1, k0); double c110 = GET_CACHE(i0+1, j0+1, k0);
  double c011 = GET_CACHE(i0, j0+1, k1); double c111 = GET_CACHE(i0+1, j0+1, k1);

  double c00 = c000 * (1.0 - rd) + c100 * rd;
  double c01 = c001 * (1.0 - rd) + c101 * rd;
  double c10 = c010 * (1.0 - rd) + c110 * rd;
  double c11 = c011 * (1.0 - rd) + c111 * rd;

  double c0 = c00 * (1.0 - td) + c10 * td;
  double c1 = c01 * (1.0 - td) + c11 * td;

  *potential_v = c0 * (1.0 - pd) + c1 * pd;
  #undef GET_CACHE
  return 0;
}
extern int cube_yy_igrf_trace(const double start[3],double target,double time,double endpoint[4]);
extern int cube_yy_footpoint_potential(double th,double ph,double *value);
static int trace_RK4(double r,double th,double ph,double target,double time,
                    double *ot,double *op,double *ratio,double *inc) {
 double start[3]={r,th,ph},out[4];int st=cube_yy_igrf_trace(start,target,time,out);
 if(!st){*ot=out[0];*op=out[1];*ratio=out[2];*inc=out[3];}return st;
}
static double interpolate_global_potential(double th,double ph) {
 double value;
 return cube_yy_footpoint_potential(th,ph,&value)?NAN:value;
}
#define GAMERA_MI_GLOBAL_SOLVE 1
static int refresh_potential_cache(double time_s) {
  const double dr=(mi_state.pot_cache_r_max-mi_state.pot_cache_r_min)/(CACHE_NR-1);
  const double dth=PI/(CACHE_NTH-1), dph=2.*PI/CACHE_NPH;
  int ranks;
  if(MPI_Comm_size(cache_comm,&ranks)!=MPI_SUCCESS)return -1;
  int owner=rank;
#ifdef SRC2_PERF_REFERENCE
  ranks=1;owner=0;
#endif
  const int count=CACHE_NR*CACHE_NTH*CACHE_NPH;
  int counts[ranks],offsets[ranks];
  for(int p=0;p<ranks;++p){offsets[p]=(int)((long long)count*p/ranks);
    counts[p]=(int)((long long)count*(p+1)/ranks)-offsets[p];}
#pragma omp parallel for schedule(dynamic,32)
  for(int flat=offsets[owner];flat<offsets[owner]+counts[owner];++flat){
    int i=flat/(CACHE_NTH*CACHE_NPH), j=(flat/CACHE_NPH)%CACHE_NTH, k=flat%CACHE_NPH;
    double th,ph,ratio,inc;
    if(trace_RK4(mi_state.pot_cache_r_min+i*dr,j*dth,k*dph,
        mi_config.ionosphere_radius_re,time_s,&th,&ph,&ratio,&inc)==0){
      mi_state.pot_cache_th_foot[flat]=th;mi_state.pot_cache_ph_foot[flat]=ph;
#ifdef GAMERA_MI_GLOBAL_SOLVE
      mi_state.pot_cache_value[flat]=interpolate_global_potential(th,ph);
#else
      mi_state.pot_cache_value[flat]=interpolate_2d_iono_potential(th,ph);
#endif
    }else{mi_state.pot_cache_th_foot[flat]=-1.;mi_state.pot_cache_value[flat]=0.;}
  }
#ifndef SRC2_PERF_REFERENCE
  double *arrays[]={mi_state.pot_cache_th_foot,mi_state.pot_cache_ph_foot,mi_state.pot_cache_value};
  for(int a=0;a<3;++a)
    if(MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DOUBLE,arrays[a],counts,offsets,
                     MPI_DOUBLE,cache_comm)!=MPI_SUCCESS)return -1;
#endif
  int unresolved = 0;
  int unresolved_by_radius[CACHE_NR] = {0};
  int unresolved_by_latitude[CACHE_NTH] = {0};
  for (int flat=0; flat<count; ++flat)
    if (mi_state.pot_cache_th_foot[flat]<0. || !isfinite(mi_state.pot_cache_value[flat])) {
      ++unresolved;
      ++unresolved_by_radius[flat/(CACHE_NTH*CACHE_NPH)];
      ++unresolved_by_latitude[(flat/CACHE_NPH)%CACHE_NTH];
    }
  if (rank==0) {
    log_info("Two caps singlefoot cache unresolved=%d count=%d radius=[%.17g,%.17g] target=%.17g", unresolved, count,
             mi_state.pot_cache_r_min, mi_state.pot_cache_r_max, mi_config.ionosphere_radius_re);
    if (unresolved)
      for (int i=0;i<CACHE_NR;++i)
        log_info("Two caps cache radial layer %d r=%.17g unresolved=%d", i,
                 mi_state.pot_cache_r_min+i*dr, unresolved_by_radius[i]);
    if (unresolved)
      for (int j=0;j<CACHE_NTH;++j)
        if (unresolved_by_latitude[j])
          log_info("Two caps cache colatitude %d theta=%.17g unresolved=%d", j,
                   j*dth, unresolved_by_latitude[j]);
  }
  if (unresolved) return -1;
#ifdef GAMERA_MI_NATIVE_WALL_MAP
  /* Retain the old volume cache solely as an explicitly inactive reference
   * during this first bounded numerical candidate. It never supplies wall E
   * or CT endpoint potentials when native-point mapping is enabled. */
  return refresh_native_wall_mapping(time_s);
#else
  return 0;
#endif
}

/* Domain comes from the actual mapped Cube wall, including angular ghosts.
 * No spherical-geometry requirement or alternative mesh is imposed. */
int cube_yy_volume_cache_configure(double low,double high,double rion) {
 if(!isfinite(low)||!isfinite(high)||!isfinite(rion)||rion<=0.||low<rion||high<=low)return -1;
 mi_state.pot_cache_r_min=low;mi_state.pot_cache_r_max=high;
 mi_config.ionosphere_radius_re=rion;mi_state.solution_ready=0;return 0;
}
void cube_yy_volume_cache_invalidate(void) { mi_state.solution_ready=0; }
int cube_yy_volume_cache_refresh(int fortran_comm,double time) {
 int local=0,global=0;
 cache_comm=MPI_Comm_f2c((MPI_Fint)fortran_comm);
 if(cache_comm==MPI_COMM_NULL)return -1;
 if(MPI_Comm_rank(cache_comm,&rank)!=MPI_SUCCESS)return -1;
 if(!isfinite(time)||time<0.||!(mi_state.pot_cache_r_max>mi_state.pot_cache_r_min))local=1;
 const size_t count=CACHE_NR*CACHE_NTH*CACHE_NPH;
 if(!mi_state.pot_cache_value)mi_state.pot_cache_value=calloc(count,sizeof(double));
 if(!mi_state.pot_cache_th_foot)mi_state.pot_cache_th_foot=calloc(count,sizeof(double));
 if(!mi_state.pot_cache_ph_foot)mi_state.pot_cache_ph_foot=calloc(count,sizeof(double));
 if(!mi_state.pot_cache_value||!mi_state.pot_cache_th_foot||!mi_state.pot_cache_ph_foot)local=1;
 if(MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MAX,cache_comm)!=MPI_SUCCESS||global)return -1;
 mi_state.solution_ready=0;
 int result=refresh_potential_cache(time);
 if(!result)mi_state.solution_ready=1;
 return result;
}
int cube_yy_volume_potential(const double point[3],double *value) {
 if(!point||!value||!isfinite(point[0])||!isfinite(point[1])||!isfinite(point[2]))return -1;
 gamera_no_vec3 p={{point[0],point[1],point[2]}};return interpolate_potential(p,value);
}
int cube_yy_volume_cache_dump(const char *path,double time) {
 FILE *f=fopen(path,"wbx");if(!f)return -1;
 int dims[3]={CACHE_NR,CACHE_NTH,CACHE_NPH};const size_t count=CACHE_NR*CACHE_NTH*CACHE_NPH;
 int bad=fwrite("YYVOL001",1,8,f)!=8 || fwrite(dims,sizeof(int),3,f)!=3 ||
  fwrite(&mi_state.pot_cache_r_min,sizeof(double),1,f)!=1 || fwrite(&mi_state.pot_cache_r_max,sizeof(double),1,f)!=1 ||
  fwrite(&time,sizeof(double),1,f)!=1 || fwrite(mi_state.pot_cache_th_foot,sizeof(double),count,f)!=count ||
  fwrite(mi_state.pot_cache_ph_foot,sizeof(double),count,f)!=count || fwrite(mi_state.pot_cache_value,sizeof(double),count,f)!=count;
 return fclose(f)||bad?-1:0;
}
