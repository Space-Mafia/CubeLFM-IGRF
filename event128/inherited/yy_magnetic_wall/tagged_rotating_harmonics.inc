static void evaluate_rotating_harmonics_direct(gamera_no_vec3 point, double time, gamera_no_vec3 *B_cart) {
    double x = point.value[0];
    double y = point.value[1];
    double z = point.value[2];
    double r = sqrt(x*x + y*y + z*z);
    
    for (int i = 0; i < GAMERA_NO_DIM; i++) B_cart->value[i] = 0.0;
    if (r < 1e-10) return;

    double theta = acos(z / r);
    double phi_grid = atan2(y, x);
    
    double omega_code = 0.0;
#ifdef GAMERA_MI_COUPLING
    omega_code = mi_config.omega_rad_s * norm_config.Time_Norm;
#endif
    double phase0 = 0.0;
#ifdef GAMERA_MI_COUPLING
    phase0 = mi_config.rotation_phase0_rad;
#endif
    double phi_eff = phi_grid - (phase0 + omega_code * time);

    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    if (sin_theta < 1e-10) sin_theta = 1e-10;

    const int degree = sh_active_degree();
    /* Geometry-only memoization. Time, rotation phase, coefficients and
     * magnetospheric fields are not cached. Exact point/norm/degree key,
     * bounded per-thread storage and uncached fallback on allocation/full. */
    typedef struct {
        gamera_no_vec3 point;
        double x_norm;
        int degree, used;
        double P[SH_N_MAX+1][SH_N_MAX+1];
        double dP[SH_N_MAX+1][SH_N_MAX+1];
        double power[SH_N_MAX+1];
    } geometry_entry;
    enum { GEOMETRY_CAP = 8192 };
    geometry_entry *geometry_table = gamera_harmonic_thread_cache(
        0, GEOMETRY_CAP * sizeof(geometry_entry));
    geometry_entry *geometry = NULL;
    int geometry_new = 1;
    uint64_t geometry_hash = UINT64_C(1469598103934665603);
    for (int a=0;a<3;++a) {
        uint64_t bits;memcpy(&bits,&point.value[a],sizeof(bits));
        geometry_hash=(geometry_hash^bits)*UINT64_C(1099511628211);
        geometry_hash^=geometry_hash>>32;
    }
    if (geometry_table) for(unsigned probe=0;probe<64;++probe) {
        geometry_entry *e=&geometry_table[(geometry_hash+probe)&(GEOMETRY_CAP-1)];
        if (!e->used) {
            e->used=1;e->point=point;e->x_norm=norm_config.x_Norm;e->degree=degree;
            geometry=e;break;
        }
        if (!memcmp(&e->point,&point,sizeof(point)) && e->x_norm==norm_config.x_Norm && e->degree==degree) {
            geometry=e;geometry_new=0;break;
        }
    }
    double local_P[SH_N_MAX+1][SH_N_MAX+1],local_dP[SH_N_MAX+1][SH_N_MAX+1];
    double local_power[SH_N_MAX+1];
    double (*P)[SH_N_MAX+1]=geometry?geometry->P:local_P;
    double (*dP)[SH_N_MAX+1]=geometry?geometry->dP:local_dP;
    double *power=geometry?geometry->power:local_power;
    if (geometry_new) {
    memset(P,0,sizeof(local_P));memset(dP,0,sizeof(local_dP));
    P[0][0] = 1.0;
    dP[0][0] = 0.0;
    
    for (int n = 1; n <= degree; n++) {
        for (int m = 0; m <= n; m++) {
            if (n == m) {
                double S = (n == 1) ? 1.0 : sqrt((double)(2 * n - 1) / (2.0 * n));
                P[n][n] = S * sin_theta * P[n-1][n-1];
                dP[n][n] = S * (cos_theta * P[n-1][n-1] + sin_theta * dP[n-1][n-1]);
            } else {
                double a = (2.0 * n - 1.0) / sqrt((double)(n * n - m * m));
                double b = sqrt((double)((n - 1) * (n - 1) - m * m) / (double)(n * n - m * m));
                P[n][m] = a * cos_theta * P[n-1][m] - b * (n > 1 ? P[n-2][m] : 0.0);
                dP[n][m] = a * (-sin_theta * P[n-1][m] + cos_theta * dP[n-1][m]) - b * (n > 1 ? dP[n-2][m] : 0.0);
            }
        }
    }

    double R_sim = (R_PLANET_KM * 1000.0) / norm_config.x_Norm;
    for(int n=1;n<=degree;++n) {
        double r_ratio=R_sim/r;
        power[n]=pow(r_ratio,n+2);
    }
    }
    double nT_to_code = 1.0e-9 / norm_config.B_Norm;
    double Br_n = 0.0, Bth_n = 0.0, Bph_n = 0.0;

#ifdef GAMERA_PERF_HARMONIC_TRIG
    /* Same argument and libm call, reused within this evaluation only.
     * No phase interpolation, angular recurrence or changed summation. */
    double cos_order[SH_N_MAX + 1], sin_order[SH_N_MAX + 1];
    unsigned char order_ready[SH_N_MAX + 1] = {0};
#endif

    for (int n = 1; n <= degree; n++) {
        double r_n2 = power[n];
        
        // ONLY sum m > 0 terms
        for (int m = 1; m <= n; m++) {
            double g = sh_coeffs_g[n][m];
            double h = sh_coeffs_h[n][m];
            if (g == 0.0 && h == 0.0) continue;

#ifdef GAMERA_PERF_HARMONIC_TRIG
            if (!order_ready[m]) {
                cos_order[m] = cos(m * phi_eff);
                sin_order[m] = sin(m * phi_eff);
                order_ready[m] = 1;
            }
            double cos_mp = cos_order[m];
            double sin_mp = sin_order[m];
#else
            double cos_mp = cos(m * phi_eff);
            double sin_mp = sin(m * phi_eff);
#endif

            Br_n += (n + 1) * r_n2 * P[n][m] * (g * cos_mp + h * sin_mp);
            Bth_n -= r_n2 * dP[n][m] * (g * cos_mp + h * sin_mp);
            Bph_n += r_n2 * P[n][m] * m * (g * sin_mp - h * cos_mp) / sin_theta;
        }
    }

    B_cart->value[0] = (Br_n * sin_theta * cos(phi_grid) + Bth_n * cos_theta * cos(phi_grid) - Bph_n * sin(phi_grid)) * nT_to_code;
    B_cart->value[1] = (Br_n * sin_theta * sin(phi_grid) + Bth_n * cos_theta * sin(phi_grid) + Bph_n * cos(phi_grid)) * nT_to_code;
    B_cart->value[2] = (Br_n * cos_theta - Bth_n * sin_theta) * nT_to_code;
}



/* Exact memoization, not phase interpolation: repeated current/old boundary
 * calls at the same fixed point reuse the original computed vector. The
 * magnetospheric perturbation is NEVER cached. Thread-local tables avoid
 * introducing races if callers later become OpenMP parallel. Bounded storage
 * falls back to direct evaluation when full; no assumption on mesh size. */
static void evaluate_rotating_harmonics(gamera_no_vec3 point, double time,
                                       gamera_no_vec3 *field) {
#ifdef SRC2_PERF_REFERENCE
    evaluate_rotating_harmonics_direct(point, time, field);
#else
    enum { CAP = 32768 };
    typedef struct {
        gamera_no_vec3 point, value[2];
        double time[2];
        unsigned used, count, next;
    } entry;
    /* Keep large storage off OpenMP thread stacks/static TLS. */
    entry *table = gamera_harmonic_thread_cache(1, CAP * sizeof(entry));
    static __thread double signature[5];
    double omega = 0.;
#ifdef GAMERA_MI_COUPLING
    omega = mi_config.omega_rad_s;
#endif
    double phase0 = 0.;
#ifdef GAMERA_MI_COUPLING
    phase0 = mi_config.rotation_phase0_rad;
#endif
    double key[5] = {norm_config.x_Norm, norm_config.B_Norm,
                     norm_config.Time_Norm, omega, phase0};
    if (!table) {
        evaluate_rotating_harmonics_direct(point, time, field);
        return;
    }
    if (memcmp(signature, key, sizeof(key))) {
        memset(table, 0, CAP * sizeof(*table));
        memcpy(signature, key, sizeof(key));
    }
    uint64_t hash = UINT64_C(1469598103934665603);
    for (int a = 0; a < 3; ++a) {
        uint64_t bits;
        memcpy(&bits, &point.value[a], sizeof(bits));
        hash = (hash ^ bits) * UINT64_C(1099511628211);
        hash ^= hash >> 32;
    }
    for (unsigned probe = 0; probe < 64; ++probe) {
        entry *e = &table[(hash + probe) & (CAP - 1)];
        if (!e->used) { e->point = point; e->used = 1; }
        if (memcmp(&e->point, &point, sizeof(point))) continue;
        for (unsigned t = 0; t < e->count; ++t)
            if (e->time[t] == time) { *field = e->value[t]; return; }
        evaluate_rotating_harmonics_direct(point, time, field);
        unsigned slot = e->next;
        e->value[slot] = *field; e->time[slot] = time;
        e->next = 1 - slot;
        if (e->count < 2) ++e->count;
        return;
    }
    evaluate_rotating_harmonics_direct(point, time, field);
#endif
}
