#include "nonorthogonal_mi_global.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int valid_config(const gamera_mi_global_solver_config *config) {
  return config != NULL && config->longitude_count >= 4U &&
         config->colatitude_count >= 5U &&
         isfinite(config->ionosphere_radius_m) &&
         config->ionosphere_radius_m > 0.0 &&
         isfinite(config->pedersen_siemens) &&
         config->pedersen_siemens > 0.0 &&
         config->maximum_iterations > 0 &&
         isfinite(config->relative_tolerance) &&
         config->relative_tolerance > 0.0 &&
         isfinite(config->absolute_tolerance) &&
         config->absolute_tolerance >= 0.0;
}

static size_t reduced_index(size_t np, size_t i, size_t j) {
  return 1U + (i - 1U) * np + j;
}

static size_t south_pole_index(size_t np, size_t nt) {
  return 1U + (nt - 2U) * np;
}

static double dot(const double *left, const double *right, size_t count) {
  double value = 0.0;
  for (size_t index = 0; index < count; ++index) {
    value += left[index] * right[index];
  }
  return value;
}

static void project_zero_sum(double *values, size_t count) {
  double mean = 0.0;
  for (size_t index = 0; index < count; ++index) {
    mean += values[index] / (double)count;
  }
  for (size_t index = 0; index < count; ++index) {
    values[index] -= mean;
  }
}

static int apply_operator(const gamera_mi_global_solver_config *config,
                          const double *potential, double *result) {
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const size_t south = south_pole_index(np, nt);
  const double pi = acos(-1.0);
  const double dtheta = pi / (double)(nt - 1U);
  const double dphi = 2.0 * pi / (double)np;
  const double pole_edge = config->pedersen_siemens * sin(0.5 * dtheta) /
                           (dtheta * dtheta);
  double north_ring_sum = 0.0;
  double south_ring_sum = 0.0;
  for (size_t j = 0; j < np; ++j) {
    north_ring_sum += potential[reduced_index(np, 1U, j)];
    south_ring_sum += potential[reduced_index(np, nt - 2U, j)];
  }
  result[0] = pole_edge * ((double)np * potential[0] - north_ring_sum);
  result[south] =
      pole_edge * ((double)np * potential[south] - south_ring_sum);
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double lower = config->pedersen_siemens *
                         sin(theta - 0.5 * dtheta) / (dtheta * dtheta);
    const double upper = config->pedersen_siemens *
                         sin(theta + 0.5 * dtheta) / (dtheta * dtheta);
    const double longitude = config->pedersen_siemens /
                             (sin(theta) * dphi * dphi);
    for (size_t j = 0; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t center = reduced_index(np, i, j);
      const double lower_value =
          i == 1U ? potential[0]
                  : potential[reduced_index(np, i - 1U, j)];
      const double upper_value =
          i + 2U == nt ? potential[south]
                       : potential[reduced_index(np, i + 1U, j)];
      result[center] =
          lower * (potential[center] - lower_value) +
          upper * (potential[center] - upper_value) +
          longitude *
              (2.0 * potential[center] -
               potential[reduced_index(np, i, jm)] -
               potential[reduced_index(np, i, jp)]);
    }
  }
  return 0;
}

static void build_diagonal(const gamera_mi_global_solver_config *config,
                           double *diagonal) {
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const size_t south = south_pole_index(np, nt);
  const double pi = acos(-1.0);
  const double dtheta = pi / (double)(nt - 1U);
  const double dphi = 2.0 * pi / (double)np;
  const double pole_edge = config->pedersen_siemens * sin(0.5 * dtheta) /
                           (dtheta * dtheta);
  diagonal[0] = (double)np * pole_edge;
  diagonal[south] = (double)np * pole_edge;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double lower = config->pedersen_siemens *
                         sin(theta - 0.5 * dtheta) / (dtheta * dtheta);
    const double upper = config->pedersen_siemens *
                         sin(theta + 0.5 * dtheta) / (dtheta * dtheta);
    const double longitude = config->pedersen_siemens /
                             (sin(theta) * dphi * dphi);
    for (size_t j = 0; j < np; ++j) {
      diagonal[reduced_index(np, i, j)] =
          lower + upper + 2.0 * longitude;
    }
  }
}

static void unpack_solution(const gamera_mi_global_solver_config *config,
                            const double *solution, double *potential) {
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const size_t south = south_pole_index(np, nt);
  double weighted_mean = 0.0;
  double weight_sum = 0.0;
  const double dtheta = acos(-1.0) / (double)(nt - 1U);
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double weight = sin((double)i * dtheta);
    for (size_t j = 0; j < np; ++j) {
      weighted_mean += weight * solution[reduced_index(np, i, j)];
      weight_sum += weight;
    }
  }
  weighted_mean = weight_sum > 0.0 ? weighted_mean / weight_sum : 0.0;
  for (size_t j = 0; j < np; ++j) {
    potential[j] = solution[0] - weighted_mean;
    potential[(nt - 1U) * np + j] = solution[south] - weighted_mean;
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0; j < np; ++j) {
      potential[i * np + j] =
          solution[reduced_index(np, i, j)] - weighted_mean;
    }
  }
}

int gamera_mi_global_solve_constant_pedersen(
    const gamera_mi_global_solver_config *config,
    const double *outward_fac_a_m2, double *potential_v,
    gamera_mi_global_solver_stats *stats) {
  if (!valid_config(config) || outward_fac_a_m2 == NULL ||
      potential_v == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  if (np > (SIZE_MAX - 2U) / (nt - 2U)) {
    return -1;
  }
  const size_t count = 2U + np * (nt - 2U);
  double *work = (double *)calloc(7U * count, sizeof(*work));
  if (work == NULL) {
    return -1;
  }
  double *right = work;
  double *diagonal = right + count;
  double *solution = diagonal + count;
  double *residual = solution + count;
  double *preconditioned = residual + count;
  double *search = preconditioned + count;
  double *operator_search = search + count;
  const double pi = acos(-1.0);
  const double dtheta = pi / (double)(nt - 1U);
  const double dphi = 2.0 * pi / (double)np;
  double fac_weighted_sum = 0.0;
  double fac_weight_sum = 0.0;
  for (size_t j = 0; j < np; ++j) {
    if (!isfinite(outward_fac_a_m2[j]) ||
        !isfinite(outward_fac_a_m2[(nt - 1U) * np + j])) {
      free(work);
      return -1;
    }
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double weight = sin((double)i * dtheta);
    for (size_t j = 0; j < np; ++j) {
      const double fac = outward_fac_a_m2[i * np + j];
      if (!isfinite(fac)) {
        free(work);
        return -1;
      }
      fac_weighted_sum += weight * fac;
      fac_weight_sum += weight;
    }
  }
  const double removed_mean = fac_weighted_sum / fac_weight_sum;
  const double area_scale = config->ionosphere_radius_m *
                            config->ionosphere_radius_m * dtheta * dphi;
  const double radius_square = config->ionosphere_radius_m *
                               config->ionosphere_radius_m;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double weight = sin((double)i * dtheta);
    for (size_t j = 0; j < np; ++j) {
      right[reduced_index(np, i, j)] =
          -radius_square * weight *
          (outward_fac_a_m2[i * np + j] - removed_mean);
    }
  }
  for (size_t j = 0; j < np; ++j) {
    const double north = potential_v[j];
    const double south = potential_v[(nt - 1U) * np + j];
    solution[0] += isfinite(north) ? north / (double)np : 0.0;
    solution[south_pole_index(np, nt)] +=
        isfinite(south) ? south / (double)np : 0.0;
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0; j < np; ++j) {
      const double value = potential_v[i * np + j];
      solution[reduced_index(np, i, j)] = isfinite(value) ? value : 0.0;
    }
  }
  project_zero_sum(solution, count);
  build_diagonal(config, diagonal);
  apply_operator(config, solution, operator_search);
  for (size_t index = 0; index < count; ++index) {
    residual[index] = right[index] - operator_search[index];
  }
  project_zero_sum(residual, count);
  const double right_norm = sqrt(fmax(0.0, dot(right, right, count)));
  const double initial = sqrt(fmax(0.0, dot(residual, residual, count)));
  const double target = fmax(config->absolute_tolerance,
                             config->relative_tolerance *
                                 fmax(right_norm, 1.0));
  int converged = initial <= target;
  int iterations = 0;
  double final = initial;
  if (!converged) {
    for (size_t index = 0; index < count; ++index) {
      preconditioned[index] = residual[index] / diagonal[index];
    }
    project_zero_sum(preconditioned, count);
    memcpy(search, preconditioned, count * sizeof(*search));
    double rho = dot(residual, preconditioned, count);
    for (int iteration = 1; iteration <= config->maximum_iterations;
         ++iteration) {
      apply_operator(config, search, operator_search);
      const double denominator = dot(search, operator_search, count);
      if (!isfinite(rho) || !isfinite(denominator) || rho <= DBL_MIN ||
          denominator <= DBL_MIN) {
        break;
      }
      const double alpha = rho / denominator;
      for (size_t index = 0; index < count; ++index) {
        solution[index] += alpha * search[index];
        residual[index] -= alpha * operator_search[index];
      }
      project_zero_sum(solution, count);
      project_zero_sum(residual, count);
      final = sqrt(fmax(0.0, dot(residual, residual, count)));
      iterations = iteration;
      if (final <= target) {
        converged = 1;
        break;
      }
      for (size_t index = 0; index < count; ++index) {
        preconditioned[index] = residual[index] / diagonal[index];
      }
      project_zero_sum(preconditioned, count);
      const double rho_new = dot(residual, preconditioned, count);
      if (!isfinite(rho_new) || rho_new <= DBL_MIN) {
        break;
      }
      const double beta = rho_new / rho;
      for (size_t index = 0; index < count; ++index) {
        search[index] = preconditioned[index] + beta * search[index];
      }
      project_zero_sum(search, count);
      rho = rho_new;
    }
  }
  if (stats != NULL) {
    stats->solver = (gamera_mi_solver_stats){
        iterations, initial, final, converged};
    stats->removed_mean_fac_a_m2 = removed_mean;
    stats->net_fac_before_a = fac_weighted_sum * area_scale;
    stats->net_fac_after_a =
        (fac_weighted_sum - removed_mean * fac_weight_sum) * area_scale;
  }
  if (!converged) {
    free(work);
    return -1;
  }
  unpack_solution(config, solution, potential_v);
  free(work);
  return 0;
}

/* ======================================================================
 * GLOBAL TENSOR SOLVER IMPLEMENTATION
 * ====================================================================== */

static int valid_global_tensor_config(const gamera_mi_global_solver_config *config,
                                      const double *pedersen,
                                      const double *hall) {
  if (config == NULL || pedersen == NULL || hall == NULL ||
      config->longitude_count < 4U || config->colatitude_count < 5U ||
      !isfinite(config->ionosphere_radius_m) ||
      !(config->ionosphere_radius_m > 0.0) ||
      config->maximum_iterations <= 0 ||
      !isfinite(config->relative_tolerance) ||
      !(config->relative_tolerance > 0.0) ||
      !isfinite(config->absolute_tolerance) ||
      config->absolute_tolerance < 0.0) {
    return 0;
  }
  if (config->longitude_count > SIZE_MAX / config->colatitude_count) return 0;
  const size_t count = config->longitude_count * config->colatitude_count;
  for (size_t index = 0U; index < count; ++index) {
    if (!isfinite(pedersen[index]) || !(pedersen[index] > 0.0) ||
        !isfinite(hall[index]) || hall[index] < 0.0) {
      return 0;
    }
  }
  return 1;
}

static double global_dipole_cos_inclination(double theta) {
  const double cosine = cos(theta);
  /* Returns negative in Northern hemisphere, positive in Southern hemisphere */
  return -2.0 * cosine / sqrt(1.0 + 3.0 * cosine * cosine);
}

static double global_inclination_square(double theta) {
  const double cosine = cos(theta);
  double inc_sq = 4.0 * cosine * cosine / (1.0 + 3.0 * cosine * cosine);
  /* Clamp to prevent singularity at the magnetic equator */
  if (inc_sq < 1.0e-6) inc_sq = 1.0e-6;
  return inc_sq;
}

static double global_tensor_f12(size_t np, const double *hall, size_t i, size_t j, double dtheta) {
  const double theta = (double)i * dtheta;
  double inclination = global_dipole_cos_inclination(theta);
  /* Clamp Hall divisor at the equator */
  if (fabs(inclination) < 1.0e-3) {
    inclination = inclination < 0.0 ? -1.0e-3 : 1.0e-3;
  }
  return -hall[i * np + j] / inclination;
}

typedef struct {
  const gamera_mi_global_solver_config *config;
  double *diagonal;
  double *lower;
  double *upper;
  double *longitude_lower;
  double *longitude_upper;
  double *pole_edge_north;
  double *pole_edge_south;
  double *preconditioner_work;
  double *storage;
} gamera_mi_global_tensor_stencil;

static void destroy_global_tensor_stencil(gamera_mi_global_tensor_stencil *stencil) {
  if (stencil != NULL) {
    free(stencil->storage);
    *stencil = (gamera_mi_global_tensor_stencil){0};
  }
}

static int build_global_tensor_stencil(
    const gamera_mi_global_solver_config *config, const double *pedersen,
    const double *hall, double *diagonal,
    gamera_mi_global_tensor_stencil *stencil) {
  if (stencil == NULL || diagonal == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const size_t count = 2U + np * (nt - 2U);
  
  double *storage = (double *)calloc(5U * count + 2U * np, sizeof(*storage));
  if (storage == NULL) {
    return -1;
  }
  
  *stencil = (gamera_mi_global_tensor_stencil){
      .config = config,
      .diagonal = diagonal,
      .lower = storage,
      .upper = storage + count,
      .longitude_lower = storage + 2U * count,
      .longitude_upper = storage + 3U * count,
      .pole_edge_north = storage + 4U * count,
      .pole_edge_south = storage + 4U * count + np,
      .preconditioner_work = storage + 4U * count + 2U * np,
      .storage = storage};
      
  const double pi = acos(-1.0);
  const double dtheta = pi / (double)(nt - 1U);
  const double dphi = 2.0 * pi / (double)np;
  const size_t south = south_pole_index(np, nt);
  
  double north_mean_pedersen = 0.0;
  double south_mean_pedersen = 0.0;
  for (size_t j = 0U; j < np; ++j) {
    north_mean_pedersen += pedersen[j] / (double)np;
    south_mean_pedersen += pedersen[(nt - 1U) * np + j] / (double)np;
  }
  
  diagonal[0] = 0.0;
  diagonal[south] = 0.0;
  
  for (size_t j = 0U; j < np; ++j) {
    double face_north = 0.5 * (north_mean_pedersen + pedersen[np + j]);
    stencil->pole_edge_north[j] = sin(0.5 * dtheta) * face_north /
         (global_inclination_square(0.5 * dtheta) * dtheta * dtheta);
    diagonal[0] += stencil->pole_edge_north[j];
    
    double face_south = 0.5 * (south_mean_pedersen + pedersen[(nt - 2U) * np + j]);
    stencil->pole_edge_south[j] = sin(pi - 0.5 * dtheta) * face_south /
         (global_inclination_square(pi - 0.5 * dtheta) * dtheta * dtheta);
    diagonal[south] += stencil->pole_edge_south[j];
  }
  
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double sine = sin(theta);
    const double lower_theta = theta - 0.5 * dtheta;
    const double upper_theta = theta + 0.5 * dtheta;
    
    for (size_t j = 0U; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t full = i * np + j;
      const size_t center = reduced_index(np, i, j);
      
      const double lower =
          sin(lower_theta) * 0.5 * (pedersen[full] + pedersen[(i - 1U) * np + j]) /
          (global_inclination_square(lower_theta) * dtheta * dtheta);
          
      const double upper =
          sin(upper_theta) * 0.5 * (pedersen[full] + pedersen[(i + 1U) * np + j]) /
          (global_inclination_square(upper_theta) * dtheta * dtheta);
          
      const double longitude_lower =
          0.5 * (pedersen[full] + pedersen[i * np + jm]) / (sine * dphi * dphi);
          
      const double longitude_upper =
          0.5 * (pedersen[full] + pedersen[i * np + jp]) / (sine * dphi * dphi);
          
      diagonal[center] = lower + upper + longitude_lower + longitude_upper;
      
      const double df12_dtheta = (global_tensor_f12(np, hall, i + 1U, j, dtheta) -
                                  global_tensor_f12(np, hall, i - 1U, j, dtheta)) / (2.0 * dtheta);
                                  
      const double df12_dphi = (global_tensor_f12(np, hall, i, jp, dtheta) -
                                global_tensor_f12(np, hall, i, jm, dtheta)) / (2.0 * dphi);
                                
      stencil->lower[center] = -lower - df12_dphi / (2.0 * dtheta);
      stencil->upper[center] = -upper + df12_dphi / (2.0 * dtheta);
      stencil->longitude_lower[center] = -longitude_lower + df12_dtheta / (2.0 * dphi);
      stencil->longitude_upper[center] = -longitude_upper - df12_dtheta / (2.0 * dphi);
    }
  }
  return 0;
}

static void apply_global_tensor_stencil(const gamera_mi_global_tensor_stencil *stencil,
                                        const double *potential_v,
                                        double *result) {
  const gamera_mi_global_solver_config *config = stencil->config;
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const size_t south = south_pole_index(np, nt);
  
  result[0] = stencil->diagonal[0] * potential_v[0];
  for (size_t j = 0U; j < np; ++j) {
    result[0] -= stencil->pole_edge_north[j] * potential_v[reduced_index(np, 1U, j)];
  }
  
  result[south] = stencil->diagonal[south] * potential_v[south];
  for (size_t j = 0U; j < np; ++j) {
    result[south] -= stencil->pole_edge_south[j] * potential_v[reduced_index(np, nt - 2U, j)];
  }
  
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0U; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t center = reduced_index(np, i, j);
      
      const double lower_value = (i == 1U) ? potential_v[0] : potential_v[reduced_index(np, i - 1U, j)];
      const double upper_value = (i + 2U == nt) ? potential_v[south] : potential_v[reduced_index(np, i + 1U, j)];
      
      result[center] =
          stencil->diagonal[center] * potential_v[center] +
          stencil->lower[center] * lower_value +
          stencil->upper[center] * upper_value +
          stencil->longitude_lower[center] * potential_v[reduced_index(np, i, jm)] +
          stencil->longitude_upper[center] * potential_v[reduced_index(np, i, jp)];
    }
  }
}

static void apply_global_tensor_preconditioner(gamera_mi_global_tensor_stencil *stencil,
                                               const double *input,
                                               double *output) {
  const gamera_mi_global_solver_config *config = stencil->config;
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const size_t south = south_pole_index(np, nt);
  double *forward = stencil->preconditioner_work;

  // Forward sweep
  forward[0] = input[0] / stencil->diagonal[0];
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0U; j < np; ++j) {
      const size_t center = reduced_index(np, i, j);
      double value = input[center];
      const double lower_val = (i == 1U) ? forward[0] : forward[reduced_index(np, i - 1U, j)];
      
      value -= stencil->lower[center] * lower_val;
      if (j > 0U) value -= stencil->longitude_lower[center] * forward[reduced_index(np, i, j - 1U)];
      if (j + 1U == np) value -= stencil->longitude_upper[center] * forward[reduced_index(np, i, 0U)];
      
      forward[center] = value / stencil->diagonal[center];
    }
  }
  
  double south_val = input[south];
  for (size_t j = 0U; j < np; ++j) {
    south_val += stencil->pole_edge_south[j] * forward[reduced_index(np, nt - 2U, j)];
  }
  forward[south] = south_val / stencil->diagonal[south];

  // Backward sweep
  output[south] = forward[south];
  for (size_t i = nt - 2U; i > 0U; --i) {
    for (size_t j_reverse = np; j_reverse > 0U; --j_reverse) {
      const size_t j = j_reverse - 1U;
      const size_t center = reduced_index(np, i, j);
      double value = stencil->diagonal[center] * forward[center];
      
      const double upper_val = (i + 2U == nt) ? output[south] : output[reduced_index(np, i + 1U, j)];
      value -= stencil->upper[center] * upper_val;
      
      if (j + 1U < np) value -= stencil->longitude_upper[center] * output[reduced_index(np, i, j + 1U)];
      if (j == 0U) value -= stencil->longitude_lower[center] * output[reduced_index(np, i, np - 1U)];
      
      output[center] = value / stencil->diagonal[center];
    }
  }
  
  double north_val = stencil->diagonal[0] * forward[0];
  for (size_t j = 0U; j < np; ++j) {
    north_val += stencil->pole_edge_north[j] * output[reduced_index(np, 1U, j)];
  }
  output[0] = north_val / stencil->diagonal[0];
}

#include "tensor_gmres.inc"
#ifdef GAMERA_MI_SPARSE_DIRECT
#include "tensor_direct.inc"
#endif

int gamera_mi_global_solve_conductance_tensor(
    const gamera_mi_global_solver_config *config, const double *fac_a_m2,
    const double *pedersen, const double *hall, double *potential_v,
    gamera_mi_global_solver_stats *stats) {
  if (stats != NULL) {
    *stats = (gamera_mi_global_solver_stats){0};
    stats->solver.initial_residual = INFINITY;
    stats->solver.final_residual = INFINITY;
  }
  if (!valid_global_tensor_config(config, pedersen, hall) || fac_a_m2 == NULL ||
      potential_v == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  if (np > (SIZE_MAX - 2U) / (nt - 2U)) {
    return -1;
  }
  const size_t count = 2U + np * (nt - 2U);
  /* Covers workspaces, including the bounded GMRES fallback. */
  if (count > SIZE_MAX / sizeof(double) / 256U) return -1;
  for (size_t k = 0; k < nt * np; ++k) {
    if (!isfinite(fac_a_m2[k]) || !isfinite(potential_v[k])) return -1;
  }
  double *work = (double *)calloc(11U * count, sizeof(double));
  if (work == NULL) return -1;
  
  double *right = work;
  double *diagonal = right + count;
  double *solution = diagonal + count;
  double *residual = solution + count;
  double *shadow = residual + count;
  double *search = shadow + count;
  double *operator_search = search + count;
  double *preconditioned_search = operator_search + count;
  double *intermediate = preconditioned_search + count;
  double *preconditioned_intermediate = intermediate + count;
  double *operator_intermediate = preconditioned_intermediate + count;

  const double pi = acos(-1.0);
  const double dtheta = pi / (double)(nt - 1U);
  const double dphi = 2.0 * pi / (double)np;
  const double radius_square = config->ionosphere_radius_m * config->ionosphere_radius_m;

  double fac_weighted_sum = 0.0;
  double fac_weight_sum = 0.0;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double weight = sin((double)i * dtheta);
    for (size_t j = 0; j < np; ++j) {
      fac_weighted_sum += weight * fac_a_m2[i * np + j];
      fac_weight_sum += weight;
    }
  }
  const double removed_mean = fac_weighted_sum / fac_weight_sum;
  
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double weight = sin((double)i * dtheta);
    for (size_t j = 0; j < np; ++j) {
      right[reduced_index(np, i, j)] = -radius_square * weight * (fac_a_m2[i * np + j] - removed_mean);
    }
  }

  const size_t south = south_pole_index(np, nt);
  for (size_t j = 0; j < np; ++j) {
    const double n_val = potential_v[j];
    const double s_val = potential_v[(nt - 1U) * np + j];
    solution[0] += isfinite(n_val) ? n_val / (double)np : 0.0;
    solution[south] += isfinite(s_val) ? s_val / (double)np : 0.0;
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0; j < np; ++j) {
      const double value = potential_v[i * np + j];
      solution[reduced_index(np, i, j)] = isfinite(value) ? value : 0.0;
    }
  }

  project_zero_sum(solution, count);

  gamera_mi_global_tensor_stencil stencil = {0};
  if (build_global_tensor_stencil(config, pedersen, hall, diagonal, &stencil) != 0) {
    free(work);
    return -1;
  }
#ifdef GAMERA_MI_SPARSE_DIRECT
  /* Alternative linear algebra only: exact same A, b and output gauge.
   * Removing the north-pole unknown fixes only the constant null mode.
   * ALL rows, including the omitted pole equation, must still pass. */
  apply_global_tensor_stencil(&stencil, solution, residual);
  for(size_t k=0;k<count;++k)residual[k]=right[k]-residual[k];
  double direct_initial=sqrt(dot(residual,residual,count));
  double direct_target=0.8*fmax(config->absolute_tolerance,
      config->relative_tolerance*fmax(sqrt(dot(right,right,count)),1.0));
  double direct_final=INFINITY;
  int direct_steps=0;
  int direct_rc=tensor_sparse_direct(&stencil,right,solution,count,
      direct_target,config->maximum_iterations,&direct_steps,&direct_final);
  int direct_ok=direct_rc==0 && isfinite(direct_final) && direct_final<=direct_target;
  if(stats){
    stats->solver=(gamera_mi_solver_stats){direct_steps,direct_initial,direct_final,direct_ok};
    stats->removed_mean_fac_a_m2=removed_mean;
    stats->net_fac_before_a=fac_weighted_sum*radius_square*dtheta*dphi;
    stats->net_fac_after_a=(fac_weighted_sum-removed_mean*fac_weight_sum)*radius_square*dtheta*dphi;
  }
  if(direct_ok)unpack_solution(config,solution,potential_v);
  destroy_global_tensor_stencil(&stencil);free(work);
  return direct_ok?0:-1;
#endif
  /* Keep Hall in A, but not in the approximate triangular inverse: the
   * antisymmetric/drift part can destroy the SGS diffusion preconditioner. */
  double *pmem = calloc(nt*np + count, sizeof(double));
  gamera_mi_global_tensor_stencil preconditioner = {0};
  if (pmem == NULL || build_global_tensor_stencil(config, pedersen, pmem,
          pmem + nt*np, &preconditioner) != 0) {
    free(pmem); destroy_global_tensor_stencil(&stencil); free(work); return -1;
  }
  
  apply_global_tensor_stencil(&stencil, solution, operator_search);
  for (size_t index = 0U; index < count; ++index) {
    residual[index] = right[index] - operator_search[index];
  }
  project_zero_sum(residual, count);
  memcpy(shadow, residual, count * sizeof(*shadow));
  
  const double right_norm = sqrt(dot(right, right, count));
  const double initial = sqrt(dot(residual, residual, count));
  /* Leave headroom for cancellation when the output gauge is unpacked and
   * independently re-applied by the caller. Never relax its requested tol. */
  const double target = 0.8 * fmax(config->absolute_tolerance, config->relative_tolerance * fmax(right_norm, 1.0));
  
  int converged = isfinite(initial) && isfinite(target) && initial <= target;
  int iterations = 0;
  double final = initial;
  double rho_previous = 1.0;
  double alpha = 1.0;
  double omega = 1.0;
  
  const int bicg_budget = config->maximum_iterations < 512 ? config->maximum_iterations : 512;
  for (int iteration = 1; !converged && iteration <= bicg_budget; ++iteration) {
    const double rho = dot(shadow, residual, count);
    if (!isfinite(rho) || fabs(rho) <= DBL_MIN || !isfinite(omega) || fabs(omega) <= DBL_MIN) break;
    
    const double beta = (rho / rho_previous) * (alpha / omega);
    for (size_t index = 0U; index < count; ++index) {
      search[index] = residual[index] + beta * (search[index] - omega * operator_search[index]);
    }
    
    apply_global_tensor_preconditioner(&stencil, search, preconditioned_search);
    project_zero_sum(preconditioned_search, count);
    apply_global_tensor_stencil(&stencil, preconditioned_search, operator_search);
    
    const double alpha_denominator = dot(shadow, operator_search, count);
    if (!isfinite(alpha_denominator) || fabs(alpha_denominator) <= DBL_MIN) break;
    
    alpha = rho / alpha_denominator;
    for (size_t index = 0U; index < count; ++index) {
      intermediate[index] = residual[index] - alpha * operator_search[index];
    }
    project_zero_sum(intermediate, count);
    
    final = sqrt(dot(intermediate, intermediate, count));
    if (final <= target) {
      for (size_t index = 0U; index < count; ++index) solution[index] += alpha * preconditioned_search[index];
      project_zero_sum(solution, count);
      iterations = iteration;
      converged = 1;
      break;
    }
    
    apply_global_tensor_preconditioner(&stencil, intermediate, preconditioned_intermediate);
    project_zero_sum(preconditioned_intermediate, count);
    apply_global_tensor_stencil(&stencil, preconditioned_intermediate, operator_intermediate);
    
    const double omega_denominator = dot(operator_intermediate, operator_intermediate, count);
    if (!isfinite(omega_denominator) || omega_denominator <= DBL_MIN) break;
    
    omega = dot(operator_intermediate, intermediate, count) / omega_denominator;
    if (!isfinite(omega) || fabs(omega) <= DBL_MIN) break;
    
    for (size_t index = 0U; index < count; ++index) {
      solution[index] += alpha * preconditioned_search[index] + omega * preconditioned_intermediate[index];
      residual[index] = intermediate[index] - omega * operator_intermediate[index];
    }
    project_zero_sum(solution, count);
    project_zero_sum(residual, count);
    
    final = sqrt(dot(residual, residual, count));
    iterations = iteration;
    
    if (final <= target || iteration % 256 == 0) {
      apply_global_tensor_stencil(&stencil, solution, operator_intermediate);
      for (size_t index = 0U; index < count; ++index) residual[index] = right[index] - operator_intermediate[index];
      project_zero_sum(residual, count);
      final = sqrt(dot(residual, residual, count));
      
      if (final <= target) {
        converged = 1;
        break;
      }
      
      memcpy(shadow, residual, count * sizeof(*shadow));
      memset(search, 0, count * sizeof(*search));
      memset(operator_search, 0, count * sizeof(*operator_search));
      rho_previous = 1.0;
      alpha = 1.0;
      omega = 1.0;
      continue;
    }
    rho_previous = rho;
  }
  
  /* Accept only the unprojected residual of the actual operator, including
   * the early BiCGSTAB exit. A recurrence residual is not a certificate. */
  apply_global_tensor_stencil(&stencil, solution, operator_intermediate);
  for (size_t k = 0; k < count; ++k)
    residual[k] = right[k] - operator_intermediate[k];
  final = sqrt(dot(residual, residual, count));
  converged = isfinite(final) && isfinite(target) && final <= target;
  if (!converged && isfinite(initial) && isfinite(target) &&
      iterations < config->maximum_iterations) {
    /* The full-tensor SGS is fast on ordinary auroral inputs but can be
     * unstable for extreme Hall gradients. Do not seed the robust fallback
     * with a diverged iterate. */
    if (!isfinite(final) || final > initial) memset(solution,0,count*sizeof(double));
    int used=tensor_gmres(&stencil,&preconditioner,right,solution,count,
                         target,config->maximum_iterations-iterations,&final);
    if(used>=0)iterations+=used;
    converged = used>=0 && isfinite(final) && final<=target;
  }
  destroy_global_tensor_stencil(&preconditioner);
  free(pmem);
  if (stats != NULL) {
    stats->solver = (gamera_mi_solver_stats){iterations, initial, final, converged};
    stats->removed_mean_fac_a_m2 = removed_mean;
    const double area_scale = radius_square * dtheta * dphi;
    stats->net_fac_before_a = fac_weighted_sum * area_scale;
    stats->net_fac_after_a = (fac_weighted_sum - removed_mean * fac_weight_sum) * area_scale;
  }
  
  if (!converged) {
    destroy_global_tensor_stencil(&stencil);
    free(work);
    return -1;
  }
  
  unpack_solution(config, solution, potential_v);
  destroy_global_tensor_stencil(&stencil);
  free(work);
  return 0;
}
