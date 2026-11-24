/* Interface adaptation only. The tagged harmonic evaluator/RK4 are unchanged. */
#define GAMERA_MI_COUPLING
#define norm_config cube_yy_norm_config
#define mi_config cube_yy_mi_config
#include "magnetic_tracer_invariants_v3.h"

Norm_t cube_yy_norm_config;
MI_config_t cube_yy_mi_config;
static int configured;

int cube_yy_igrf_configure(double length_m, double magnetic_tesla,
                          double omega_rad_s, double phase_rad) {
    if (!isfinite(length_m) || length_m <= 0.0 ||
        !isfinite(magnetic_tesla) || magnetic_tesla <= 0.0 ||
        !isfinite(omega_rad_s) || !isfinite(phase_rad)) return -1;
    norm_config.x_Norm = length_m;
    norm_config.B_Norm = magnetic_tesla;
    mi_config.omega_rad_s = omega_rad_s;
    mi_config.rotation_phase0_rad = phase_rad;
    cube_yy_prepare_invariants();
    cube_yy_field_cache_reset();
    configured = 1;
    return 0;
}

int cube_yy_igrf_field(const double xyz[3], double time_s,
                      double direction[3], double *magnitude) {
    if (!configured || !xyz || !direction || !magnitude || !isfinite(time_s)) return -1;
    if (!isfinite(xyz[0]) || !isfinite(xyz[1]) || !isfinite(xyz[2]) ||
        (xyz[0] == 0.0 && xyz[1] == 0.0 && xyz[2] == 0.0)) return -1;
    get_B_direction_cartesian(xyz[0],xyz[1],xyz[2],time_s,1.0,
        &direction[0],&direction[1],&direction[2],magnitude);
    return isfinite(direction[0]) && isfinite(direction[1]) &&
           isfinite(direction[2]) && isfinite(*magnitude) && *magnitude > 0.0 ? 0 : -1;
}

int cube_yy_igrf_trace(const double spherical_start[3], double target_radius,
                      double time_s, double endpoint[4]) {
    if (!configured || !spherical_start || !endpoint) return -1;
    /* The caller must use this status; an escaped/failed trace is not a zero FAC. */
    return trace_RK4(spherical_start[0],spherical_start[1],spherical_start[2],
        target_radius,time_s,&endpoint[0],&endpoint[1],&endpoint[2],&endpoint[3]);
}
