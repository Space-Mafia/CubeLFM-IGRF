/* Thin CT entrypoints. Both numerical bodies below are exact tagged text. */
#include <math.h>
#include <stddef.h>
#define GAMERA_MI_COUPLING
#define norm_config cube_yy_ct_norm_config
#define mi_config cube_yy_ct_mi_config
#include "config.h"
#include "background_coeffs.h"
#define GAMERA_NO_DIM 3
typedef struct { double value[3]; } gamera_no_vec3;
Norm_t cube_yy_ct_norm_config;
MI_config_t cube_yy_ct_mi_config;
static int configured;

#include "tagged_harmonic_body.inc"
#include "tagged_cubic_taper.inc"

int cube_yy_ct_configure(double length_m, double magnetic_tesla,
                        double time_unit_s, double omega_rad_s, double phase_rad) {
    if (!isfinite(length_m) || length_m <= 0.0 ||
        !isfinite(magnetic_tesla) || magnetic_tesla <= 0.0 ||
        !isfinite(time_unit_s) || time_unit_s <= 0.0 ||
        !isfinite(omega_rad_s) || !isfinite(phase_rad)) return -1;
    norm_config.x_Norm=length_m;
    norm_config.B_Norm=magnetic_tesla;
    norm_config.Time_Norm=time_unit_s;
    mi_config.omega_rad_s=omega_rad_s;
    mi_config.rotation_phase0_rad=phase_rad;
    configured=1;
    return 0;
}

int cube_yy_ct_evaluate(const double xyz[3], double time_code, int m0_only,
                       int tagged_taper, double field[3], double potential[3]) {
    if (!configured || !xyz || !field || !potential || !isfinite(time_code) ||
        (m0_only != 0 && m0_only != 1) ||
        (tagged_taper != 0 && tagged_taper != 1)) return -1;
    for (int i=0;i<3;++i) if (!isfinite(xyz[i])) return -1;
    double radius=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
    if (!(radius>0.0)) return -1;
    gamera_no_vec3 point={{xyz[0],xyz[1],xyz[2]}},b,a;
    evaluate_magnetic_harmonics(point,time_code,m0_only,&b,&a);
    /* Exact released16..24 model-length taper, separately from uncut MI field. */
    double w=tagged_taper ? cubic_ramp_down(radius,16.0,8.0) : 1.0;
    for (int i=0;i<3;++i) {
        field[i]=w*b.value[i];potential[i]=w*a.value[i];
        if (!isfinite(field[i]) || !isfinite(potential[i])) return -1;
    }
    return 0;
}
