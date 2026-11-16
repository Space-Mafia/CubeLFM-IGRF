/* Exact released four-point interpolation arithmetic. Only the coordinate
 * provider and sparse MPI accessor are replaced by saved Cube logical data.
 * Radial interpolation composes two unchanged angular interpolations. This
 * is a diagnostic; no source mask, stencil retune or new Poisson/tracer. */
#include <math.h>
static const double *saved;
static struct { int nj_global,nk_global; } config;
static int global_patch_value(int layer,int face,int j,int k,double *v) {
 if(!saved||!v||layer<0||layer>=6||face<0||face>=6||j<0||j>=config.nj_global||k<0||k>=config.nk_global)return -1;
 *v=saved[layer+6*(j+config.nj_global*(k+config.nk_global*face))];return isfinite(*v)?0:-1;
}
static void interpolation_coordinate(double logical, int count, int *base,
                                     double *fraction) {
  if (logical <= 0.0) {
    *base = 0;
    *fraction = 0.0;
  } else if (logical >= (double)(count - 1)) {
    *base = count - 2;
    *fraction = 1.0;
  } else {
    *base = (int)floor(logical);
    *fraction = logical - (double)*base;
  }
}
static int interpolate_patch_source(int source, int patch,
                                    const double point[3], double *value) {
  if(!point||!value)return -1;
  const double q_theta=point[1],q_phi=point[2];
  int base_theta, base_phi;
  double fraction_theta, fraction_phi;
  interpolation_coordinate(q_theta, config.nj_global, &base_theta,
                           &fraction_theta);
  interpolation_coordinate(q_phi, config.nk_global, &base_phi,
                           &fraction_phi);
  *value = 0.0;
  for (int dj = 0; dj <= 1; ++dj) {
    for (int dk = 0; dk <= 1; ++dk) {
      double corner;
      if (global_patch_value(source, patch, base_theta + dj, base_phi + dk,
                             &corner) != 0) {
        return -1;
      }
      const double wj = dj == 0 ? 1.0 - fraction_theta : fraction_theta;
      const double wk = dk == 0 ? 1.0 - fraction_phi : fraction_phi;
      *value += wj * wk * corner;
    }
  }
  return isfinite(*value) ? 0 : -1;
}

int cube_released_ghost_scalar_interpolate(const double logical[3],int face,int nang,
 const double *parallel,double *value,int *clamped) {
 if(!logical||!parallel||!value||!clamped)return -1;
 for(int d=0;d<3;++d)if(!isfinite(logical[d]))return -1;
 config.nj_global=nang+2;config.nk_global=nang+2;
 saved=parallel;int radial;double fraction,lo,hi;
 interpolation_coordinate(logical[0],6,&radial,&fraction);
 *clamped=(logical[0]<0.||logical[0]>5.?1:0)|
           (logical[1]<0.||logical[1]>(double)(config.nj_global-1)?2:0)|
           (logical[2]<0.||logical[2]>(double)(config.nk_global-1)?4:0);
 if(interpolate_patch_source(radial,face,logical,&lo)||
    interpolate_patch_source(radial+1,face,logical,&hi))return -1;
 *value=(1.-fraction)*lo+fraction*hi;return isfinite(*value)?0:-1;
}
