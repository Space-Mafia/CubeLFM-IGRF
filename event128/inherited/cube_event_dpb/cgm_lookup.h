#ifdef GAMERA_EVENT_IGRF
#include <event_cgm_table.h>
#else
#include "cgm_table.h"
#endif
static void cgm_interp(const double *a,int nr,int nc,double u,double v,double q[3]) {
  u=fmax(0.,fmin(nr-1.,u));int i=(int)floor(u);if(i==nr-1)i--;
  double w=u-i;v=fmod(v,nc);if(v<0)v+=nc;int j=(int)floor(v),jp=(j+1)%nc;double z=v-j;
  double n=0.;for(int k=0;k<3;k++){
    q[k]=(1-w)*((1-z)*a[3*(i*nc+j)+k]+z*a[3*(i*nc+jp)+k])+
            w*((1-z)*a[3*((i+1)*nc+j)+k]+z*a[3*((i+1)*nc+jp)+k]);n+=q[k]*q[k];
  }
  n=sqrt(n);for(int k=0;k<3;k++)q[k]/=n;
}
static double cgm_sun_longitude(double time_s) {
  double phase=mi_rotation_phase_rad(time_s),c=cos(phase),s=sin(phase);
  double sb[3]={c*mi_config.sunward_direction[0]+s*mi_config.sunward_direction[1],
             -s*mi_config.sunward_direction[0]+c*mi_config.sunward_direction[1],mi_config.sunward_direction[2]};
  double a=0.,b=0.;for(int k=0;k<3;k++){a+=sb[k]*cgm_Q[k*3];b+=sb[k]*cgm_Q[k*3+1];}
  return atan2(b,a);
}
static void cgm_foot(int h,double lat,double phi_mlt,double time_s,double alpha,double q[3]) {
  double v[3];
  cgm_interp(cgm_feet+h*CGM_NLAT*CGM_NLON*3,CGM_NLAT,CGM_NLON,
             (lat-50.)/.5,(phi_mlt-PI+alpha)*180./PI,v);
  double phase=mi_rotation_phase_rad(time_s),c=cos(phase),s=sin(phase);
  q[0]=c*v[0]-s*v[1];q[1]=s*v[0]+c*v[1];q[2]=v[2];
}
static double cgm_mask(double th,double ph,double time_s,double alpha,size_t np) {
  double v[3];
  cgm_interp(cgm_inverse,CGM_NTH,CGM_NLON,th*180./PI,
            (ph-mi_rotation_phase_rad(time_s))*180./PI,v);
  double lat=asin(fmin(1.,fabs(v[2])))*180./PI;
  if(lat<50.)return 0.;if(lat>82.)return 1.;
  int h=v[2]>=0?0:1;
  double j=fmod((atan2(v[1],v[0])-alpha+PI)/(2*PI)*np,(double)np);if(j<0)j+=np;
  size_t lo=(size_t)floor(j),hi=(lo+1)%np;double w=j-lo;
  double edge=(1-w)*mi_state.dpb_boundary_latitude_deg[h][lo]+w*mi_state.dpb_boundary_latitude_deg[h][hi];
  double u=fmax(0.,fmin(1.,(lat-edge)/mi_config.dpb_transition_width_deg));
  return u*u*(3.-2.*u);
}
