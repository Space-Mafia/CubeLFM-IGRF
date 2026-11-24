#ifndef GAMERA_AURORAL_DPB_MASK_H
#define GAMERA_AURORAL_DPB_MASK_H
#include <math.h>

/* Apply once to freshly computed precipitation, irrespective of its diffuse
 * or accelerated classification. Equal flux factors preserve mean energy;
 * the existing Robinson P/H laws therefore scale with sqrt(mask).
 * EUV, FAC, conductance floors and time filtering are handled by the caller. */
static int gamera_apply_auroral_dpb_mask(double mask, double *number_flux,
    double *energy_flux, double *auroral_p, double *auroral_h) {
  if (!isfinite(mask)) return -1;
  const double m = fmax(0., fmin(1., mask));
  const double s = sqrt(m);
  *number_flux *= m;
  *energy_flux *= m;
  *auroral_p *= s;
  *auroral_h *= s;
  return 0;
}
#endif
