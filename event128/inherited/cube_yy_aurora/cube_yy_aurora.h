#ifndef CUBE_YY_AURORA_H
#define CUBE_YY_AURORA_H

typedef struct {
    double euv_p, euv_h, auroral_p, auroral_h;
    double number_flux_cm2_s, energy_flux_erg_cm2_s, mean_energy_kev;
    double pedersen, hall;
    int diffuse;
} cube_auroral_point;

int cube_yy_aurora_point(double time_s, double zenith_rad, int trace_valid,
                         double density_g_cm3, double sound_cm_s,
                         double parallel_iono_a_m2, double signed_inclination_cos,
                         double dpb_mask, int history_ready, double previous_p,
                         double previous_h, cube_auroral_point *result);
#endif
