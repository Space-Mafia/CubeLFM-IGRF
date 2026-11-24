#ifndef MAGNETIC_TRACER_H
#define MAGNETIC_TRACER_H

#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "common.h"
#include "config.h"
#include "background_coeffs.h" // Swapped from JRM33

/* Algebraically invariant recurrence factors. Recomputed by the existing
 * explicit configure call, with the original expressions and FP64 order.
 * Coefficients, active degree and normalization are still the tagged inputs.
 */
static int cube_yy_cached_degree;
static double cube_yy_s[SH_N_MAX+1];
static double cube_yy_a[SH_N_MAX+1][SH_N_MAX+1],cube_yy_b[SH_N_MAX+1][SH_N_MAX+1];
static double cube_yy_g[SH_N_MAX+1][SH_N_MAX+1],cube_yy_h[SH_N_MAX+1][SH_N_MAX+1];
static void cube_yy_prepare_invariants(void) {
    cube_yy_cached_degree=sh_active_degree();
    for(int n=1;n<=cube_yy_cached_degree;++n){
        cube_yy_s[n]=(n==1)?1.0:sqrt((double)(2*n-1)/(2.0*n));
        for(int m=0;m<=n;++m){
            if(m<n){
                cube_yy_a[n][m]=(2.0*n-1.0)/sqrt((double)(n*n-m*m));
                cube_yy_b[n][m]=sqrt((double)((n-1)*(n-1)-m*m)/(double)(n*n-m*m));
            }
            cube_yy_g[n][m]=sh_coeffs_g[n][m]*1e-9/norm_config.B_Norm;
            cube_yy_h[n][m]=sh_coeffs_h[n][m]*1e-9/norm_config.B_Norm;
        }
    }
}

// -------------------------------------------------------------------------
// 1. DYNAMIC FIELD EVALUATOR (WITH COROTATION)
// -------------------------------------------------------------------------
static inline void get_planetary_B(double r_sim, double th, double ph_grid, double time_s,
                                   double *Br, double *Bth, double *Bph) {
    *Br = 0.0; *Bth = 0.0; *Bph = 0.0;

    // --- APPLY COROTATION ---
   // --- APPLY COROTATION ---
    double rotation_phase = 0.0;
#ifdef GAMERA_MI_COUPLING
    rotation_phase = mi_rotation_phase_rad(time_s);
#endif

    double ph_eff = ph_grid - rotation_phase;

    double cos_th = cos(th);
    double sin_th = sin(th);
    if (fabs(sin_th) < 1e-15) sin_th = (sin_th < 0) ? -1e-15 : 1e-15;

    // Robust Schmidt semi-normalized Legendre polynomials
    const int degree = cube_yy_cached_degree;
    double P[degree + 1][degree + 1];
    double dP[degree + 1][degree + 1];
    memset(P,0,sizeof(P));memset(dP,0,sizeof(dP));
    P[0][0] = 1.0;
    dP[0][0] = 0.0;

    for (int n = 1; n <= degree; n++) {
        for (int m = 0; m <= n; m++) {
            if (n == m) {
                double S = cube_yy_s[n];
                P[n][n] = S * sin_th * P[n-1][n-1];
                dP[n][n] = S * (cos_th * P[n-1][n-1] + sin_th * dP[n-1][n-1]);
            } else {
                double a = cube_yy_a[n][m];
                double b = cube_yy_b[n][m];
                P[n][m] = a * cos_th * P[n-1][m] - b * (n > 1 ? P[n-2][m] : 0.0);
                dP[n][m] = a * (-sin_th * P[n-1][m] + cos_th * dP[n-1][m]) - b * (n > 1 ? dP[n-2][m] : 0.0);
            }
        }
    }

    double r_ratio = (R_PLANET_KM * 1000.0) / (r_sim * norm_config.x_Norm);

    // Precompute radial decay
    double r_decay[SH_N_MAX + 3];
    r_decay[0] = 1.0;
    for (int n = 1; n <= degree + 2; n++) {
        r_decay[n] = r_decay[n-1] * r_ratio;
    }

    for (int m = 0; m <= degree; m++) {
        double cos_mph = cos(m * ph_eff);
        double sin_mph = sin(m * ph_eff);
        double m_over_sin_th = (m > 0) ? (m / sin_th) : 0.0;

        for (int n = (m == 0 ? 1 : m); n <= degree; n++) {
            double g = cube_yy_g[n][m];
            double h = cube_yy_h[n][m];

            if (g == 0.0 && h == 0.0) continue;

            double r_n2 = r_decay[n + 2];

            *Br += (n + 1) * r_n2 * P[n][m] * (g * cos_mph + h * sin_mph);
            *Bth -= r_n2 * dP[n][m] * (g * cos_mph + h * sin_mph);
            if (m > 0) {
                *Bph += r_n2 * m_over_sin_th * P[n][m] * (g * sin_mph - h * cos_mph);
            }
        }
    }
}

// -------------------------------------------------------------------------
// 2. CARTESIAN WRAPPER FOR TRACING
// -------------------------------------------------------------------------
static inline void get_B_direction_cartesian_uncached(double x, double y, double z, double time_s, double trace_dir,
                                             double *Bx, double *By, double *Bz, double *B_mag) {
    double r = sqrt(x*x + y*y + z*z);
    double th = acos(z / r);
    double ph = atan2(y, x);
    if (ph < 0) ph += 2.0 * PI;

    double Br, Bth, Bph;
    get_planetary_B(r, th, ph, time_s, &Br, &Bth, &Bph);

    double sth = sin(th), cth = cos(th);
    double sph = sin(ph), cph = cos(ph);

    // Reconstruct Cartesian using the TRUE SM grid longitude (ph)
    double Bx_raw = Br * sth * cph + Bth * cth * cph - Bph * sph;
    double By_raw = Br * sth * sph + Bth * cth * sph + Bph * cph;
    double Bz_raw = Br * cth       - Bth * sth;

    *B_mag = sqrt(Bx_raw*Bx_raw + By_raw*By_raw + Bz_raw*Bz_raw);

    *Bx = trace_dir * (Bx_raw / *B_mag);
    *By = trace_dir * (By_raw / *B_mag);
    *Bz = trace_dir * (Bz_raw / *B_mag);
}

// -------------------------------------------------------------------------
// 3. RK4 TRACER
// -------------------------------------------------------------------------

/* Exact bounded memoization. Eviction only recomputes the unchanged
 * evaluator; it cannot change any field, tolerance, status or trace step.
 * In particular shell-bisection traces share identical RK4 prefix queries.
 */
#define CUBE_YY_FIELD_CACHE_SIZE 16384u
struct cube_yy_field_entry { uint64_t key[4],generation; double value[4]; };
static struct cube_yy_field_entry *cube_yy_field_cache;
static uint64_t cube_yy_field_generation=1,cube_yy_field_time_word;
static int cube_yy_field_time_valid,cube_yy_field_allocation_attempted;
static void cube_yy_field_cache_reset(void) {
    cube_yy_field_time_valid=0;
    if(++cube_yy_field_generation==0){
        if(cube_yy_field_cache)memset(cube_yy_field_cache,0,
            CUBE_YY_FIELD_CACHE_SIZE*sizeof(*cube_yy_field_cache));
        cube_yy_field_generation=1;
    }
}
static inline void get_B_direction_cartesian(double x,double y,double z,
    double time_s,double trace_dir,double *bx,double *by,double *bz,double *mag) {
    uint64_t tw,key[4];
    if(!isfinite(x)||!isfinite(y)||!isfinite(z)||!isfinite(time_s)||!isfinite(trace_dir)){
        get_B_direction_cartesian_uncached(x,y,z,time_s,trace_dir,bx,by,bz,mag);return;
    }
    if(!cube_yy_field_allocation_attempted){
        cube_yy_field_cache=calloc(CUBE_YY_FIELD_CACHE_SIZE,sizeof(*cube_yy_field_cache));
        cube_yy_field_allocation_attempted=1;
    }
    if(!cube_yy_field_cache){
        get_B_direction_cartesian_uncached(x,y,z,time_s,trace_dir,bx,by,bz,mag);return;
    }
    memcpy(&tw,&time_s,sizeof(tw));
    if(!cube_yy_field_time_valid||tw!=cube_yy_field_time_word){
        cube_yy_field_cache_reset();cube_yy_field_time_word=tw;cube_yy_field_time_valid=1;
    }
    memcpy(key,&x,8);memcpy(key+1,&y,8);memcpy(key+2,&z,8);memcpy(key+3,&trace_dir,8);
    uint64_t h=UINT64_C(1469598103934665603);
    for(int i=0;i<4;++i){h=(h^key[i])*UINT64_C(1099511628211);h^=h>>32;}
    struct cube_yy_field_entry *e=&cube_yy_field_cache[h&(CUBE_YY_FIELD_CACHE_SIZE-1u)];
    if(e->generation==cube_yy_field_generation && memcmp(e->key,key,sizeof(key))==0){
        *bx=e->value[0];*by=e->value[1];*bz=e->value[2];*mag=e->value[3];return;
    }
    get_B_direction_cartesian_uncached(x,y,z,time_s,trace_dir,bx,by,bz,mag);
    if(isfinite(*bx)&&isfinite(*by)&&isfinite(*bz)&&isfinite(*mag)){
        memcpy(e->key,key,sizeof(key));e->value[0]=*bx;e->value[1]=*by;e->value[2]=*bz;e->value[3]=*mag;
        e->generation=cube_yy_field_generation;
    }
}

/* Refine only the final RK4 step: evaluate the same field-line ODE from
 * the step's starting point, not a straight-line projection onto the shell. */
static inline void gamera_trace_substep(double x, double y, double z,
    double h, double time_s, double direction, double *qx, double *qy, double *qz) {
    double ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,bmag;
    get_B_direction_cartesian(x,y,z,time_s,direction,&ax,&ay,&az,&bmag);
    get_B_direction_cartesian(x+.5*h*ax,y+.5*h*ay,z+.5*h*az,time_s,direction,&bx,&by,&bz,&bmag);
    get_B_direction_cartesian(x+.5*h*bx,y+.5*h*by,z+.5*h*bz,time_s,direction,&cx,&cy,&cz,&bmag);
    get_B_direction_cartesian(x+h*cx,y+h*cy,z+h*cz,time_s,direction,&dx,&dy,&dz,&bmag);
    *qx=x+h*(ax+2*bx+2*cx+dx)/6.;
    *qy=y+h*(ay+2*by+2*cy+dy)/6.;
    *qz=z+h*(az+2*bz+2*cz+dz)/6.;
}

static inline int trace_RK4(double r_start, double th_start, double ph_start, double r_target, double time_s,
                            double *th_foot, double *ph_foot, double *B_ratio, double *cos_inclination) {

    if (!th_foot || !ph_foot || !B_ratio || !cos_inclination ||
        !isfinite(r_start) || !isfinite(r_target) || !isfinite(th_start) ||
        !isfinite(ph_start) || !isfinite(time_s) || !(r_start > 0.) ||
        !(r_target > 0.) || th_start < 0. || th_start > PI) return -1;
    double x = r_start * sin(th_start) * cos(ph_start);
    double y = r_start * sin(th_start) * sin(ph_start);
    double z = r_start * cos(th_start);

    double r_inv = 1.0 / r_start;
    double nx = x * r_inv, ny = y * r_inv, nz = z * r_inv;

    double Bx_start, By_start, Bz_start, B_mag_start;
    get_B_direction_cartesian(x, y, z, time_s, 1.0, &Bx_start, &By_start, &Bz_start, &B_mag_start);

    if (!(B_mag_start > 0.) || !isfinite(B_mag_start)) return -1;
    double B_radial = Bx_start * nx + By_start * ny + Bz_start * nz;
    if (r_target == r_start) {
        *th_foot=th_start; *ph_foot=fmod(ph_start,2.*PI);
        if (*ph_foot < 0.) *ph_foot += 2.*PI;
        *B_ratio=1.; *cos_inclination=fabs(B_radial); return 0;
    }
    double trace_dir = 1.0;
    if (r_target < r_start) {
        trace_dir = (B_radial <= 0.0) ? 1.0 : -1.0;
    } else {
        trace_dir = (B_radial >= 0.0) ? 1.0 : -1.0;
    }

    double ds = (r_target < r_start) ? 0.04 : 0.08;

    int max_steps = 40000;
    int step = 0;
    double r_current = r_start;

    while (step < max_steps) {
        const double x_before=x,y_before=y,z_before=z,r_before=r_current;
        double k1x, k1y, k1z, bmag1;
        get_B_direction_cartesian(x, y, z, time_s, trace_dir, &k1x, &k1y, &k1z, &bmag1);

        double k2x, k2y, k2z, bmag2;
        get_B_direction_cartesian(x + 0.5*ds*k1x, y + 0.5*ds*k1y, z + 0.5*ds*k1z, time_s, trace_dir, &k2x, &k2y, &k2z, &bmag2);

        double k3x, k3y, k3z, bmag3;
        get_B_direction_cartesian(x + 0.5*ds*k2x, y + 0.5*ds*k2y, z + 0.5*ds*k2z, time_s, trace_dir, &k3x, &k3y, &k3z, &bmag3);

        double k4x, k4y, k4z, bmag4;
        get_B_direction_cartesian(x + ds*k3x, y + ds*k3y, z + ds*k3z, time_s, trace_dir, &k4x, &k4y, &k4z, &bmag4);

        x += (ds / 6.0) * (k1x + 2.0*k2x + 2.0*k3x + k4x);
        y += (ds / 6.0) * (k1y + 2.0*k2y + 2.0*k3y + k4y);
        z += (ds / 6.0) * (k1z + 2.0*k2z + 2.0*k3z + k4z);

        r_current = sqrt(x*x + y*y + z*z);
        if (r_current < 0.95) return -1;

        if ((r_target < r_start && r_current <= r_target) ||
            (r_target > r_start && r_current >= r_target)) {

            /* Land on the requested shell. The old code used an overshot
             * point, introducing step-count sawteeth into footpoint mapping. */
            double lo=0.,hi=ds,flo=r_before-r_target,fhi=r_current-r_target;
            int landed=0;
            for (int refine=0;refine<40;++refine) {
                double h=(lo*fhi-hi*flo)/(fhi-flo);
                if (!(h>lo && h<hi) || refine%8==7) h=.5*(lo+hi);
                gamera_trace_substep(x_before,y_before,z_before,h,time_s,trace_dir,&x,&y,&z);
                r_current=sqrt(x*x+y*y+z*z);
                const double f=r_current-r_target;
                if (!isfinite(f)) return -1;
                if (fabs(f)<=1.e-12*fmax(1.,r_target)) { landed=1; break; }
                if ((f>0.)==(flo>0.)) {lo=h;flo=f;} else {hi=h;fhi=f;}
            }
            if (!landed) return -1;
            *th_foot = acos(fmax(-1.,fmin(1.,z / r_current)));
            *ph_foot = atan2(y, x);
            if (*ph_foot < 0.0) *ph_foot += 2.0 * PI;

            /* Ratio is evaluated below at the actual target shell. */

            double Bx_foot, By_foot, Bz_foot, Bmag_foot;
            get_B_direction_cartesian(x, y, z, time_s, 1.0, &Bx_foot, &By_foot, &Bz_foot, &Bmag_foot);
            if (!(Bmag_foot>0.) || !isfinite(Bmag_foot)) return -1;
            *B_ratio = B_mag_start / Bmag_foot;
            double r_inv_foot = 1.0 / r_current;
            double radial_dot = Bx_foot * (x * r_inv_foot) + By_foot * (y * r_inv_foot) + Bz_foot * (z * r_inv_foot);
            *cos_inclination = fabs(radial_dot);

            return 0;
        }
        step++;
    }
    return -1;
}
#endif // MAGNETIC_TRACER_H
