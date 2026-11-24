// background_coeffs.h (Test 1: 10-Degree Tilted Dipole toward +X)
#ifndef BACKGROUND_COEFFS_H
#define BACKGROUND_COEFFS_H

#define SH_N_MAX 14
#define R_PLANET_KM 6371.2 

#ifdef GAMERA_EVENT_IGRF
#include <event_igrf_coeffs.h>
#else
static const double sh_coeffs_g[SH_N_MAX + 1][SH_N_MAX + 1] = {
    {0.0},                                                            // n=0
    {-29404.8, -1450.9},                                              // n=1: g_1^0, g_1^1
    {-2499.6, 2982.0, 1677.0},                                        // n=2: g_2^0, g_2^1, g_2^2
    {1363.2, -2381.2, 1236.2, 525.7},                                 // n=3: g_3^0 ... g_3^3
    {903.0, 809.5, 86.3, -309.4, 48.0},                               // n=4: g_4^0 ... g_4^4
    {-234.3, 363.2, 187.8, -140.7, -151.2, 13.5},                     // n=5: g_5^0 ... g_5^5
    {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}     // n=6 through 14
};

static const double sh_coeffs_h[SH_N_MAX + 1][SH_N_MAX + 1] = {
    {0.0},                                                            // n=0
    {0.0, 4652.5},                                                    // n=1: h_1^0, h_1^1
    {0.0, -2991.6, -734.6},                                           // n=2: h_2^0, h_2^1, h_2^2
    {0.0, -82.1, 241.9, -543.4},                                      // n=3: h_3^0 ... h_3^3
    {0.0, 282.1, -158.4, 199.7, -349.7},                              // n=4: h_4^0 ... h_4^4
    {0.0, 47.7, 208.3, -121.2, 32.3, 99.1},                           // n=5: h_5^0 ... h_5^5
    {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}, {0.0}     // n=6 through 14
};
#endif

/* All nonzero coefficients are retained; trailing exactly-zero degrees need
 * no Legendre recurrence. This also adapts if the coefficient table changes. */
static inline int sh_active_degree(void) {
#ifdef SRC2_PERF_REFERENCE
    return SH_N_MAX;
#else
    for (int n = SH_N_MAX; n > 0; --n)
        for (int m = 0; m <= n; ++m)
            if (sh_coeffs_g[n][m] != 0.0 || sh_coeffs_h[n][m] != 0.0) return n;
    return 0;
#endif
}
#endif // BACKGROUND_COEFFS_H
