#include "solar_wind.h"
#include <math.h>
#include <stddef.h>

/* One immutable event input per MPI rank. Donor coordinates put the Sun
   toward -X, while the Cube nose is +X. A pi rotation about the common
   spin axis transforms both polar and axial vectors with determinant +1. */
static gamera_solar_wind_series wind;
static int ready;

int cube_event_wind_open(const char *path)
{
    const gamera_solar_wind_hdf5_units si = {
        .physical_units = 1, .time_norm = 1, .density_norm = 1,
        .velocity_norm = 1, .pressure_norm = 1, .magnetic_norm = 1,
        .velocity_si_scale = 1000
    };
    if (ready || !path) return -1;
    gamera_solar_wind_init(&wind);
    if (gamera_solar_wind_load_hdf5(path, &si, &wind)) return -1;
    if (wind.enforce_bx_relation || wind.count < 2) {
        gamera_solar_wind_destroy(&wind);
        return -1;
    }
    wind.front_along_velocity = 1;
    wind.linear_interpolation = 1;
    ready = 1;
    return 0;
}

int cube_event_wind_sample(double t, const double x_cube_cm[3],
                           double gamma, double q[8])
{
    gamera_solar_wind_state sample;
    double delay;
    if (!ready || !x_cube_cm || !q || !isfinite(gamma) || gamma <= 1 ||
        !isfinite(t) || t < wind.time[0] || t > wind.time[wind.count-1]) return -1;
    const double p[3] = {-0.01*x_cube_cm[0], -0.01*x_cube_cm[1], 0.01*x_cube_cm[2]};
    if (gamera_solar_wind_sample_at(&wind, p, t, &sample, &delay)) return -1;
    const double epoch = t - delay;
    if (!isfinite(epoch) || epoch < wind.time[0] || epoch > wind.time[wind.count-1] ||
        sample.density <= 0 || sample.pressure <= 0) return -1;
    q[0] = sample.density * 1.e-3;
    q[4] = 100.0 * sqrt(gamma*sample.pressure/sample.density);
    for (int i=0;i<3;i++) {
        double sign = i==2 ? 1.0 : -1.0;
        q[i+1] = sign * 100.0 * sample.velocity[i];
        q[i+5] = sign * 1.e4 * sample.magnetic[i];
    }
    for (int i=0;i<8;i++) if (!isfinite(q[i])) return -1;
    return 0;
}
