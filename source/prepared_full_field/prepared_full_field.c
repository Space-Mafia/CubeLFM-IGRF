/* Fixed-position spatial work reused across changing rotation epochs.
 * Exact donor expressions/order, all degrees and coefficients unchanged.
 * Internal donor RK4 calls remain in the original compiled object. */
#define GAMERA_MI_COUPLING
#define norm_config cube_yy_norm_config
#define mi_config cube_yy_mi_config
#include "magnetic_tracer_invariants_v3.h"
extern Norm_t cube_yy_norm_config;
extern MI_config_t cube_yy_mi_config;
extern int __real_cube_yy_igrf_configure(double,double,double,double);
extern int __real_cube_yy_igrf_field(const double*,double,double*,double*);
enum { PREPARED_CAP=32768 };
typedef struct {
 uint64_t key[3],time_word;
 int used,value_ready;
 double phi,sth,cth,sph,cph,sin_poly;
 double direction[3],magnitude;
 double term_br[15],term_bt[15],term_bp[15],term_g[15],term_h[15],axis_br,axis_bt;
 int term_m[15],term_count;
 double *data;
} spatial_entry;
static spatial_entry *table;
static double config_signature[4];
static int configured,allocation_tried;
static uint64_t hits,misses,evaluations,fallbacks;
static void clear_geometry(void) {
 if(table)for(unsigned i=0;i<PREPARED_CAP;++i){free(table[i].data);table[i].data=NULL;table[i].used=0;table[i].value_ready=0;}
 hits=misses=evaluations=fallbacks=0;
}
int __wrap_cube_yy_igrf_configure(double l,double b,double w,double p) {
 int st=__real_cube_yy_igrf_configure(l,b,w,p);
 if(st)return st;
 clear_geometry();cube_yy_prepare_invariants();
 config_signature[0]=l;config_signature[1]=b;config_signature[2]=w;config_signature[3]=p;configured=1;
 return 0;
}
static int prepare(spatial_entry *e,const double xyz[3]) {
 const double x=xyz[0],y=xyz[1],z=xyz[2];
 const double r=sqrt(x*x+y*y+z*z),th=acos(z/r);
 double ph=atan2(y,x);if(ph<0)ph+=2.*PI;
 e->phi=ph;e->sth=sin(th);e->cth=cos(th);e->sph=sin(ph);e->cph=cos(ph);
 e->sin_poly=e->sth;if(fabs(e->sin_poly)<1.e-15)e->sin_poly=(e->sin_poly<0)?-1.e-15:1.e-15;
 const int degree=cube_yy_cached_degree,width=degree+1;
 const size_t count=(size_t)2*width*width+SH_N_MAX+3;
 double *data=calloc(count,sizeof(double));if(!data)return 1;
 double *P=data,*dP=P+width*width,*r_decay=dP+width*width;
 P[0]=1.;dP[0]=0.;
 for(int n=1;n<=degree;++n)for(int m=0;m<=n;++m){
  const size_t q=n*width+m,old=(n-1)*width+m;
  if(n==m){
   const double S=cube_yy_s[n];
   P[q]=S*e->sin_poly*P[(n-1)*width+n-1];
   dP[q]=S*(e->cth*P[(n-1)*width+n-1]+e->sin_poly*dP[(n-1)*width+n-1]);
  }else{
   const double a=cube_yy_a[n][m],b=cube_yy_b[n][m];
   P[q]=a*e->cth*P[old]-b*(n>1?P[(n-2)*width+m]:0.);
   dP[q]=a*(-e->sin_poly*P[old]+e->cth*dP[old])-b*(n>1?dP[(n-2)*width+m]:0.);
  }
 }
 const double r_ratio=(R_PLANET_KM*1000.)/(r*norm_config.x_Norm);
 r_decay[0]=1.;for(int n=1;n<=degree+2;++n)r_decay[n]=r_decay[n-1]*r_ratio;
 e->axis_br=0.;e->axis_bt=0.;e->term_count=0;
 if(degree==5){
  for(int m=0;m<=degree;++m){
   const double ms=(m>0)?(m/e->sin_poly):0.;
   for(int n=(m==0?1:m);n<=degree;++n){
    const double g=cube_yy_g[n][m],h=cube_yy_h[n][m];if(g==0.&&h==0.)continue;
    const double rn=r_decay[n+2];
    if(m==0){
     e->axis_br+=(n+1)*rn*P[n*width+m]*(g*1.+h*0.);
     e->axis_bt-=rn*dP[n*width+m]*(g*1.+h*0.);
    }else{
     const int t=e->term_count++;e->term_m[t]=m;e->term_g[t]=g;e->term_h[t]=h;
     e->term_br[t]=(n+1)*rn*P[n*width+m];
     e->term_bt[t]=rn*dP[n*width+m];
     e->term_bp[t]=rn*ms*P[n*width+m];
    }
   }
  }
 }
 e->data=data;return 0;
}
static void evaluate(spatial_entry *e,double time_s,double direction[3],double *magnitude) {
 const int degree=cube_yy_cached_degree,width=degree+1;
 const double *P=e->data,*dP=P+width*width,*r_decay=dP+width*width;
 const double rotation_phase=mi_rotation_phase_rad(time_s),ph_eff=e->phi-rotation_phase;
 double Br=0.,Bth=0.,Bph=0.;
 if(degree==5){
  Br=e->axis_br;Bth=e->axis_bt;
  double cm[6],sm[6];for(int m=1;m<=degree;++m){cm[m]=cos(m*ph_eff);sm[m]=sin(m*ph_eff);}
  for(int t=0;t<e->term_count;++t){
   const int m=e->term_m[t];const double g=e->term_g[t],h=e->term_h[t];
   const double common=g*cm[m]+h*sm[m];
   Br+=e->term_br[t]*common;Bth-=e->term_bt[t]*common;
   Bph+=e->term_bp[t]*(g*sm[m]-h*cm[m]);
  }
 }else{
 for(int m=0;m<=degree;++m){
  const double cos_mph=cos(m*ph_eff),sin_mph=sin(m*ph_eff);
  const double m_over_sin_th=(m>0)?(m/e->sin_poly):0.;
  for(int n=(m==0?1:m);n<=degree;++n){
   const double g=cube_yy_g[n][m],h=cube_yy_h[n][m];
   if(g==0.&&h==0.)continue;
   const double r_n2=r_decay[n+2];
   Br+=(n+1)*r_n2*P[n*width+m]*(g*cos_mph+h*sin_mph);
   Bth-=r_n2*dP[n*width+m]*(g*cos_mph+h*sin_mph);
   if(m>0)Bph+=r_n2*m_over_sin_th*P[n*width+m]*(g*sin_mph-h*cos_mph);
  }
 }
 }
 const double bx=Br*e->sth*e->cph+Bth*e->cth*e->cph-Bph*e->sph;
 const double by=Br*e->sth*e->sph+Bth*e->cth*e->sph+Bph*e->cph;
 const double bz=Br*e->cth-Bth*e->sth;
 *magnitude=sqrt(bx*bx+by*by+bz*bz);
 direction[0]=1.*(bx/ *magnitude);direction[1]=1.*(by/ *magnitude);direction[2]=1.*(bz/ *magnitude);
}
int __wrap_cube_yy_igrf_field(const double xyz[3],double time_s,double direction[3],double *magnitude) {
 const double sig[4]={norm_config.x_Norm,norm_config.B_Norm,mi_config.omega_rad_s,mi_config.rotation_phase0_rad};
 if(!configured||!xyz||!direction||!magnitude||!isfinite(time_s)||
    !isfinite(xyz[0])||!isfinite(xyz[1])||!isfinite(xyz[2])||
    (xyz[0]==0.&&xyz[1]==0.&&xyz[2]==0.)||memcmp(sig,config_signature,sizeof(sig)))
  return __real_cube_yy_igrf_field(xyz,time_s,direction,magnitude);
 if(!allocation_tried){allocation_tried=1;table=calloc(PREPARED_CAP,sizeof(*table));}
 if(!table){++fallbacks;return __real_cube_yy_igrf_field(xyz,time_s,direction,magnitude);}
 uint64_t key[3],tw,h=UINT64_C(1469598103934665603);memcpy(key,xyz,sizeof(key));memcpy(&tw,&time_s,8);
 for(int i=0;i<3;++i){h=(h^key[i])*UINT64_C(1099511628211);h^=h>>32;}
 spatial_entry *e=NULL;
 for(unsigned probe=0;probe<64;++probe){
  spatial_entry *candidate=&table[(h+probe)&(PREPARED_CAP-1)];
  if(!candidate->used){
   if(prepare(candidate,xyz))break;
   memcpy(candidate->key,key,sizeof(key));candidate->used=1;e=candidate;++misses;break;
  }
  if(!memcmp(candidate->key,key,sizeof(key))){e=candidate;++hits;break;}
 }
 if(!e){++fallbacks;return __real_cube_yy_igrf_field(xyz,time_s,direction,magnitude);}
 if(e->value_ready&&e->time_word==tw){memcpy(direction,e->direction,24);*magnitude=e->magnitude;return 0;}
 evaluate(e,time_s,direction,magnitude);++evaluations;
 if(!isfinite(direction[0])||!isfinite(direction[1])||!isfinite(direction[2])||!isfinite(*magnitude)||!(*magnitude>0.))
  return __real_cube_yy_igrf_field(xyz,time_s,direction,magnitude);
 memcpy(e->direction,direction,24);e->magnitude=*magnitude;e->time_word=tw;e->value_ready=1;
 return 0;
}
void cube_yy_prepared_field_metrics(uint64_t out[4]) {out[0]=hits;out[1]=misses;out[2]=evaluations;out[3]=fallbacks;}
