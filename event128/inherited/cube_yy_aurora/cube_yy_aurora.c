#include "cube_yy_aurora.h"
#include "nonorthogonal_electron_precipitation.h"
#include "auroral_dpb_mask.h"
#include <math.h>

/* This event selects the coefficients from the released two-cap case.
   parallel_iono is j_parallel(source)*B_iono/B_source in A/m^2,
   before multiplication by inclination to form the radial MI source. */
int cube_yy_aurora_point(double time_s, double zenith_rad, int trace_valid,
                         double density_g_cm3, double sound_cm_s,
                         double parallel_iono_a_m2, double signed_inclination_cos,
                         double dpb_mask, int history_ready, double previous_p,
                         double previous_h, cube_auroral_point *result)
{
    if (!result || !isfinite(time_s) || time_s < 0.0 ||
        !isfinite(dpb_mask) || dpb_mask < 0.0 || dpb_mask > 1.0 ||
        (trace_valid != 0 && trace_valid != 1) ||
        (history_ready != 0 && history_ready != 1)) return -1;
    cube_auroral_point point = {0};
    if (gamera_mi_lompe_euv_conductance(zenith_rad, 150.0, &point.euv_p,
                                        &point.euv_h)) return -1;
    if (trace_valid) {
        if (!isfinite(signed_inclination_cos)) return -1;
        gamera_mi_fedder95_config config = gamera_mi_fedder95_default_config();
        config.beta = 0.8724646;
        config.hall_model = GAMERA_MI_HALL_KAEPPLER_2015;
        gamera_mi_electron_precipitation precip;
        int polarity = signed_inclination_cos < 0.0 ? -1 : 1;
        if (gamera_mi_fedder95_electron_precipitation(
                &config, density_g_cm3 * 1000.0, sound_cm_s * 0.01,
                parallel_iono_a_m2, fmax(1.0, point.euv_p),
                fmin(1.0, time_s/600.0), polarity, &precip)) return -1;
        point.number_flux_cm2_s = precip.number_flux_cm2_s;
        point.energy_flux_erg_cm2_s = precip.energy_flux_erg_cm2_s;
        point.mean_energy_kev = precip.average_energy_kev;
        point.auroral_p = precip.pedersen_siemens;
        point.auroral_h = precip.hall_siemens;
        point.diffuse = precip.is_diffuse;
    }
    if (gamera_apply_auroral_dpb_mask(dpb_mask, &point.number_flux_cm2_s,
            &point.energy_flux_erg_cm2_s, &point.auroral_p, &point.auroral_h)) return -1;
    if (gamera_mi_combine_conductance(point.euv_p, point.euv_h,
            point.auroral_p, point.auroral_h, 1.0, 0.0, 6.0,
            &point.pedersen, &point.hall)) return -1;
    if (history_ready) {
        if (!isfinite(previous_p) || !isfinite(previous_h) ||
            previous_p <= 0.0 || previous_h < 0.0) return -1;
        const double alpha = 1.0 - exp(-10.0 / 30.0);
        point.pedersen = previous_p + alpha * (point.pedersen - previous_p);
        point.hall = previous_h + alpha * (point.hall - previous_h);
    }
    *result = point;
    return 0;
}
