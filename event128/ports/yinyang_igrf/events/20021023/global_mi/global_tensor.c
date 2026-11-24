#include "nonorthogonal_mi_global.h"
#include <math.h>

int gamera_mi_global_solve_tensor_checked(const gamera_mi_global_solver_config *,
    const double *, const double *, const double *, double *, gamera_mi_global_solver_stats *);

/* One closed-sphere solve, with the historical YYG gauge and all-row check. */
int cube_event_global_solve(int np, int nt, double radius, const double *fac,
    const double *pedersen, const double *hall, double *phi, double *metrics, int *steps)
{
    const gamera_mi_global_solver_config cfg = {(size_t)np, (size_t)nt, radius,
        5., 20000, 1.e-10, 1.e-8};
    gamera_mi_global_solver_stats s = {0};
    if (gamera_mi_global_solve_tensor_checked(&cfg,fac,pedersen,hall,phi,&s))
        return 1;
    metrics[0] = s.solver.final_residual;
    metrics[1] = s.removed_mean_fac_a_m2;
    metrics[2] = s.net_fac_before_a;
    metrics[3] = s.net_fac_after_a;
    *steps = s.solver.iterations;
    return 0;
}
