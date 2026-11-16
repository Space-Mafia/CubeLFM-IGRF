static void evaluate_magnetic_harmonics(gamera_no_vec3 point, double time, int m0_only, 
                                        gamera_no_vec3 *B_cart, gamera_no_vec3 *A_cart) {
    double x = point.value[0];
    double y = point.value[1];
    double z = point.value[2];
    double r = sqrt(x*x + y*y + z*z);
    
    // Initialize outputs
    for (int i = 0; i < GAMERA_NO_DIM; i++) { 
        B_cart->value[i] = 0.0; 
        if (A_cart != NULL) A_cart->value[i] = 0.0; 
    }
    if (r < 1e-10) return;

    double theta = acos(z / r);
    double phi_grid = atan2(y, x);
    
    // --- APPLY COROTATION ---
    // Extract rotation rate from MI Config
    double omega_code = 0.0;
#ifdef GAMERA_MI_COUPLING
    omega_code = mi_config.omega_rad_s * norm_config.Time_Norm;
#endif
    
    // Rotate the longitude into the geographic frame for evaluation
    double phase0 = 0.0;
#ifdef GAMERA_MI_COUPLING
    phase0 = mi_config.rotation_phase0_rad;
#endif
    double phi_eff = phi_grid - (phase0 + omega_code * time);

    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    if (sin_theta < 1e-10) sin_theta = 1e-10;

    // Calculate Schmidt semi-normalized Legendre polynomials and derivatives
    const int degree = sh_active_degree();
    double P[SH_N_MAX + 1][SH_N_MAX + 1] = {0.0};
    double dP[SH_N_MAX + 1][SH_N_MAX + 1] = {0.0};
    P[0][0] = 1.0;
    dP[0][0] = 0.0;
    
    for (int n = 1; n <= degree; n++) {
        for (int m = 0; m <= n; m++) {
            if (n == m) {
                // Diagonal terms
                double S = (n == 1) ? 1.0 : sqrt((double)(2 * n - 1) / (2.0 * n));
                P[n][n] = S * sin_theta * P[n-1][n-1];
                dP[n][n] = S * (cos_theta * P[n-1][n-1] + sin_theta * dP[n-1][n-1]);
            } else {
                // Off-diagonal terms
                double a = (2.0 * n - 1.0) / sqrt((double)(n * n - m * m));
                double b = sqrt((double)((n - 1) * (n - 1) - m * m) / (double)(n * n - m * m));
                
                P[n][m] = a * cos_theta * P[n-1][m] - b * (n > 1 ? P[n-2][m] : 0.0);
                dP[n][m] = a * (-sin_theta * P[n-1][m] + cos_theta * dP[n-1][m]) - b * (n > 1 ? dP[n-2][m] : 0.0);
            }
        }
    }

    double nT_to_code = 1.0e-9 / norm_config.B_Norm;
    double R_sim = (R_PLANET_KM * 1000.0) / norm_config.x_Norm;

    for (int n = 1; n <= degree; n++) {
        double r_ratio = R_sim / r;
        double r_n2 = pow(r_ratio, n + 2);
        
        double Br_n = 0.0, Bth_n = 0.0, Bph_n = 0.0;

        for (int m = 0; m <= n; m++) {
            if (m0_only && m > 0) continue;
            
            double g = sh_coeffs_g[n][m];
            double h = sh_coeffs_h[n][m];
            
            if (g == 0.0 && h == 0.0) continue;

            double cos_mp = cos(m * phi_eff);
            double sin_mp = sin(m * phi_eff);

            Br_n += (n + 1) * r_n2 * P[n][m] * (g * cos_mp + h * sin_mp);
            Bth_n -= r_n2 * dP[n][m] * (g * cos_mp + h * sin_mp);
            if (m > 0) {
                Bph_n += r_n2 * P[n][m] * m * (g * sin_mp - h * cos_mp) / sin_theta;
            }
        }

        // Convert degree-n Spherical B to Cartesian B using the TRUE physical grid longitude (phi_grid)
        double Bx_n = (Br_n * sin_theta * cos(phi_grid) + Bth_n * cos_theta * cos(phi_grid) - Bph_n * sin(phi_grid)) * nT_to_code;
        double By_n = (Br_n * sin_theta * sin(phi_grid) + Bth_n * cos_theta * sin(phi_grid) + Bph_n * cos(phi_grid)) * nT_to_code;
        double Bz_n = (Br_n * cos_theta - Bth_n * sin_theta) * nT_to_code;

        B_cart->value[0] += Bx_n;
        B_cart->value[1] += By_n;
        B_cart->value[2] += Bz_n;

        // Cartesian Vector Potential (calculated correctly per degree n)
        if (A_cart != NULL) {
             A_cart->value[0] += (y * Bz_n - z * By_n) / (double)n;
             A_cart->value[1] += (z * Bx_n - x * Bz_n) / (double)n;
             A_cart->value[2] += (x * By_n - y * Bx_n) / (double)n;
        }
    }
}
