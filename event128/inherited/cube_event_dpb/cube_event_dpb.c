#include <math.h>
#include <stddef.h>
#include "event_model_config.h"
#include "nonorthogonal_mi_dpb.h"

#define PI 3.14159265358979323846
/* Geographic frame maps donor x,y into Cube by a proper rotation of pi. */
static const struct {
    double sunward_direction[3], dpb_transition_width_deg;
} mi_config = {{-GAMERA_EVENT_SUN_X,-GAMERA_EVENT_SUN_Y,GAMERA_EVENT_SUN_Z},3.0};
static struct { const double *dpb_boundary_latitude_deg[2]; } mi_state;

static double mi_rotation_phase_rad(double t)
{
    return GAMERA_EVENT_PHASE0_RAD + PI + 0.000072921150*t;
}
#include "cgm_lookup.h"

void cube_event_dpb_feet(int nt, int np, double time, double *points)
{
    const double alpha=cgm_sun_longitude(time);
    for(int h=0;h<2;h++) for(int i=0;i<nt;i++) for(int j=0;j<np;j++){
        const double latitude=90.-90.*i/(nt-1);
        cgm_foot(h,latitude,2.*PI*j/np,time,alpha,points+3*((h*nt+i)*np+j));
    }
}

int cube_event_dpb_step(int nt,int np,const double *fac,const double *pot,
                        int *ready,double *radius,double *edge,double *mask)
{
    const gamera_mi_hybrid_dpb_config config=gamera_mi_hybrid_dpb_default_config();
    for(int h=0;h<2;h++){
        gamera_mi_hybrid_dpb_stats stats={0};
        gamera_mi_hybrid_dpb_state state={ready[h],radius[h]};
        if(gamera_mi_hybrid_dpb_update(&config,nt,np,PI/2.,40.*PI/180.,
               fac+h*nt*np,1.,pot+h*nt*np,10.,&state,edge+h*np,mask+h*nt*np,&stats))return -1;
        radius[h]=state.radius_deg;ready[h]=state.initialized;
    }
    return 0;
}

int cube_event_geographic_mask(int n,int np,double time,const double *feet,
                               const double *edge,double *mask,double *zenith,double *mlat,double *mlt)
{
    mi_state.dpb_boundary_latitude_deg[0]=edge;
    mi_state.dpb_boundary_latitude_deg[1]=edge+np;
    const double alpha=cgm_sun_longitude(time);
    for(int k=0;k<n;k++){
        const double *v=feet+3*k;
        double z=0;
        for(int i=0;i<3;i++)z+=v[i]*mi_config.sunward_direction[i];
        zenith[k]=acos(fmin(1.,fmax(-1.,z)));
        const double th=atan2(hypot(v[0],v[1]),v[2]),ph=atan2(v[1],v[0]);
        double magnetic[3];
        cgm_interp(cgm_inverse,CGM_NTH,CGM_NLON,th*180./PI,(ph-mi_rotation_phase_rad(time))*180./PI,magnetic);
        mlat[k]=asin(fmax(-1.,fmin(1.,magnetic[2])))*180./PI;
        mlt[k]=fmod(36.+(atan2(magnetic[1],magnetic[0])-alpha)*12./PI,24.);
        mask[k]=cgm_mask(th,ph,time,alpha,np);
        if(!isfinite(mask[k])||!isfinite(z)||!isfinite(th)||!isfinite(ph))return -1;
    }
    return 0;
}
