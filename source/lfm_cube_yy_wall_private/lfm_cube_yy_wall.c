/* Namespace/units host for unchanged released YY wall routines below. */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#define GAMERA_NO_DIM 3
typedef struct { double value[3]; } gamera_no_vec3;
static struct { double x_Norm,u_Norm,B_Norm,Time_Norm; } norm_config;
static struct { double omega_rad_s; } mi_config;
static struct { int solution_ready; } mi_state;
static double time_sim;
/* Derived data only. The released electric/drift arithmetic remains below. */
static gamera_no_vec3 cube_yy_observed_background,cube_yy_observed_drift;
struct cube_yy_wall_entry {
 uint64_t key[3],generation,time_word;
 gamera_no_vec3 electric,background,drift;
 int drift_valid;
};
#define CUBE_YY_WALL_CACHE_SIZE 4096u
static struct cube_yy_wall_entry cube_yy_wall_cache[CUBE_YY_WALL_CACHE_SIZE];
static uint64_t cube_yy_wall_generation=1;
static int cube_yy_wall_cache_enabled=1;
static int64_t cube_yy_wall_counts[4];
void cube_yy_wall_cache_metrics(int64_t counts[4]){memcpy(counts,cube_yy_wall_counts,sizeof(cube_yy_wall_counts));}
void cube_yy_wall_invalidate(void){
 memset(cube_yy_wall_counts,0,sizeof(cube_yy_wall_counts));
 if(++cube_yy_wall_generation==0){memset(cube_yy_wall_cache,0,sizeof(cube_yy_wall_cache));cube_yy_wall_generation=1;}
}
void cube_yy_wall_cache_enable(int enabled){cube_yy_wall_cache_enabled=enabled!=0;cube_yy_wall_invalidate();}

extern int cube_yy_igrf_field(const double xyz[3],double time_s,double direction[3],double *magnitude);
extern int cube_yy_potential_callback(const double xyz[3],double *potential_v);
static int interpolate_potential(gamera_no_vec3 point,double *phi) {
 return cube_yy_potential_callback(point.value,phi);
}
static void get_B_direction_cartesian(double x,double y,double z,double time_s,double sign,
 double *bx,double *by,double *bz,double *mag) {
 double xyz[3]={x,y,z},direction[3];
 if(cube_yy_igrf_field(xyz,time_s,direction,mag)) {
  *bx=*by=*bz=*mag=NAN;return;
 }
 *bx=sign*direction[0];*by=sign*direction[1];*bz=sign*direction[2];
}
static double vector_dot(gamera_no_vec3 left, gamera_no_vec3 right) {
  double result = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result += left.value[component] * right.value[component];
  }
  return result;
}

static gamera_no_vec3 vector_cross(gamera_no_vec3 left,
                                   gamera_no_vec3 right) {
  return (gamera_no_vec3){{
      left.value[1] * right.value[2] - left.value[2] * right.value[1],
      left.value[2] * right.value[0] - left.value[0] * right.value[2],
      left.value[0] * right.value[1] - left.value[1] * right.value[0]}};
}

static double vector_norm(gamera_no_vec3 value) {
  return sqrt(fmax(0.0, vector_dot(value, value)));
}

int gamera_mi_compose_ghost_velocity(gamera_no_vec3 source_velocity,
                                     gamera_no_vec3 background_magnetic,
                                     gamera_no_vec3 drift_velocity,
                                     gamera_no_vec3 *ghost_velocity) {
  if (ghost_velocity == NULL) {
    return -1;
  }
  const double background_square =
      vector_dot(background_magnetic, background_magnetic);
  if (!(background_square > 0.0) || !isfinite(background_square)) {
    return -1;
  }
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    ghost_velocity->value[component] =
        2.0 * drift_velocity.value[component] -
        source_velocity.value[component];
    if (!isfinite(ghost_velocity->value[component])) {
      return -1;
    }
  }
  return 0;
}

static int volume_cached_wall_electric(gamera_no_vec3 point,
                                        gamera_no_vec3 *electric) {
  if (electric == NULL || !mi_state.solution_ready) {
    return -1;
  }
  const double radius = vector_norm(point);
  if (!(radius > 0.0)) {
    return -1;
  }
  const double epsilon = 1.0e-4 * fmax(1.0, radius);
  const double electric_norm = norm_config.u_Norm * norm_config.B_Norm;
  
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    gamera_no_vec3 plus = point;
    gamera_no_vec3 minus = point;
    plus.value[component] += epsilon;
    minus.value[component] -= epsilon;
    double plus_potential, minus_potential;
    if (interpolate_potential(plus, &plus_potential) != 0 ||
        interpolate_potential(minus, &minus_potential) != 0) {
      return -1;
    }
    electric->value[component] =
        -(plus_potential - minus_potential) /
        (2.0 * epsilon * norm_config.x_Norm * electric_norm);
  }
  return 0;
}

static int compose_wall_from_electric(gamera_no_vec3 point,
                                       gamera_no_vec3 source_velocity,
                                       gamera_no_vec3 electric,
                                       gamera_no_vec3 *ghost_velocity) {
  const double radius=vector_norm(point);
  if (!ghost_velocity || !mi_state.solution_ready || !(radius>0.)) return -1;
  gamera_no_vec3 background;
  
  // Use the tracer's full evaluator to get the true, fully-rotated planetary field
  double time_s = time_sim * norm_config.Time_Norm;
  double Bx, By, Bz, Bmag;
  get_B_direction_cartesian(point.value[0], point.value[1], point.value[2], time_s, 1.0, &Bx, &By, &Bz, &Bmag);
  
  background.value[0] = Bx * Bmag;
  background.value[1] = By * Bmag;
  background.value[2] = Bz * Bmag;

  const double background_square = vector_dot(background, background);
  if (!(background_square > 0.0)) {
    return -1;
  }
  
  gamera_no_vec3 drift = vector_cross(electric, background);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    drift.value[component] /= background_square;
  }
  
  // --- INJECT COROTATION ---
  const double omega_code = mi_config.omega_rad_s * norm_config.Time_Norm;
  drift.value[0] -= omega_code * point.value[1];
  drift.value[1] += omega_code * point.value[0];
  // -------------------------

  const gamera_no_vec3 radial = {{point.value[0] / radius,
                                  point.value[1] / radius,
                                  point.value[2] / radius}};
  const double drift_radial = vector_dot(drift, radial);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    drift.value[component] -= drift_radial * radial.value[component];
  }
  
  cube_yy_observed_background=background;cube_yy_observed_drift=drift;
  return gamera_mi_compose_ghost_velocity(source_velocity, background, drift,
                                          ghost_velocity);
}

int cube_yy_wall_configure(double model_length_m,double omega_rad_s) {
 if(!(model_length_m>0.)||!isfinite(model_length_m)||!isfinite(omega_rad_s))return -1;
 norm_config.x_Norm=model_length_m;norm_config.u_Norm=1.;norm_config.B_Norm=1.e-9;
 norm_config.Time_Norm=model_length_m;mi_config.omega_rad_s=omega_rad_s;mi_state.solution_ready=1;
 cube_yy_wall_invalidate();
 return 0;
}
int cube_yy_wall_ghost(const double position_m[3],double time_s,
 const double source_m_s[3],double ghost_m_s[3]) {
 gamera_no_vec3 point,source,ghost,electric;
 time_sim=time_s/norm_config.Time_Norm;
 for(int i=0;i<3;++i){point.value[i]=position_m[i]/norm_config.x_Norm;source.value[i]=source_m_s[i];}
 if(!cube_yy_wall_cache_enabled){
  if(volume_cached_wall_electric(point,&electric))return -1;
  if(compose_wall_from_electric(point,source,electric,&ghost))return -1;
 }else{
  uint64_t key[3],tw,h=UINT64_C(1469598103934665603);
  memcpy(key,point.value,24);memcpy(&tw,&time_s,8);
  for(int i=0;i<3;++i){h=(h^key[i])*UINT64_C(1099511628211);h^=h>>32;}
  struct cube_yy_wall_entry *entry=&cube_yy_wall_cache[h&(CUBE_YY_WALL_CACHE_SIZE-1u)];
  if(entry->generation!=cube_yy_wall_generation || memcmp(entry->key,key,24)){
   ++cube_yy_wall_counts[1];
   if(volume_cached_wall_electric(point,&electric))return -1;
   memcpy(entry->key,key,24);entry->generation=cube_yy_wall_generation;
   entry->electric=electric;entry->drift_valid=0;
  }else {++cube_yy_wall_counts[0];electric=entry->electric;}
  if(entry->drift_valid && entry->time_word==tw){
   ++cube_yy_wall_counts[2];
   if(gamera_mi_compose_ghost_velocity(source,entry->background,entry->drift,&ghost))return -1;
  }else{
   ++cube_yy_wall_counts[3];
   if(compose_wall_from_electric(point,source,electric,&ghost))return -1;
   entry->time_word=tw;entry->background=cube_yy_observed_background;entry->drift=cube_yy_observed_drift;
   entry->drift_valid=1;
  }
 }
 for(int i=0;i<3;++i)ghost_m_s[i]=ghost.value[i];
 return 0;
}
