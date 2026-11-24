#ifndef GAMERA_NONORTHOGONAL_MI_GLOBAL_H
#define GAMERA_NONORTHOGONAL_MI_GLOBAL_H

#include "nonorthogonal_mi.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  size_t longitude_count;
  size_t colatitude_count; /* full sphere, including both poles */
  double ionosphere_radius_m;
  double pedersen_siemens;
  int maximum_iterations;
  double relative_tolerance;
  double absolute_tolerance;
} gamera_mi_global_solver_config;

typedef struct {
  gamera_mi_solver_stats solver;
  /* Area-weighted outward FAC mean removed to satisfy the closed-sphere
   * compatibility condition.  A positive value is outward from Earth. */
  double removed_mean_fac_a_m2;
  double net_fac_before_a;
  double net_fac_after_a;
} gamera_mi_global_solver_stats;

/*
 * Solve -div_Omega(Sigma_P grad_Omega Phi) = -R^2 j_out on a closed sphere.
 *
 * The longitude coordinate is periodic, both geographic poles are single
 * degrees of freedom, and no artificial equatorial boundary is introduced.
 * Since the operator has one constant null mode, the solver removes the
 * area-weighted mean of j_out and returns the zero-area-mean potential gauge.
 * Arrays are theta-major, longitude-contiguous, theta in [0,pi].  Pole values
 * are duplicated across longitude on input/output for convenient I/O.
 */
int gamera_mi_global_solve_constant_pedersen(
    const gamera_mi_global_solver_config *config,
    const double *outward_fac_a_m2, double *potential_v,
    gamera_mi_global_solver_stats *stats);

int gamera_mi_global_solve_conductance_tensor(
    const gamera_mi_global_solver_config *config, 
    const double *outward_fac_a_m2,
    const double *pedersen_s, 
    const double *hall_s, 
    double *potential_v,
    gamera_mi_global_solver_stats *stats);

#ifdef __cplusplus
}
#endif

#endif
