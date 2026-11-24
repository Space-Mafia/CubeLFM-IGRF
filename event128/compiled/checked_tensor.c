/* Independent residual check of the repaired solver; supplied physical
 * tensor coefficients remain unchanged. */
#include "nonorthogonal_mi_global.c"
#include <stdio.h>
int gamera_mi_global_solve_tensor_checked(
    const gamera_mi_global_solver_config *c, const double *fac,
    const double *p, const double *h, double *phi,
    gamera_mi_global_solver_stats *stats) {
  if(!valid_global_tensor_config(c,p,h) || !fac || !phi || !stats)return -1;
  size_t nt=c->colatitude_count,np=c->longitude_count,n=2+np*(nt-2);
  for(size_t k=0;k<nt*np;++k)
    if(!isfinite(fac[k]) || !isfinite(p[k]) || !isfinite(h[k]) || !isfinite(phi[k])) return -1;
  int rc=gamera_mi_global_solve_conductance_tensor(c,fac,p,h,phi,stats);
  if(rc || !stats->solver.converged || !isfinite(stats->solver.final_residual))return -1;
  for(size_t k=0;k<nt*np;++k)if(!isfinite(phi[k]))return -1;
  double *mem=calloc(3*n,sizeof(double));if(!mem)return -1;
  double *diag=mem,*x=mem+n,*res=mem+2*n;
  gamera_mi_global_tensor_stencil s={0};
  if(build_global_tensor_stencil(c,p,h,diag,&s)){free(mem);return -1;}
  x[0]=phi[0];x[n-1]=phi[(nt-1)*np];
  for(size_t i=1;i+1<nt;++i)for(size_t j=0;j<np;++j)x[reduced_index(np,i,j)]=phi[i*np+j];
  apply_global_tensor_stencil(&s,x,res);
  double rn=0.,bn=0.;
  for(size_t i=1;i+1<nt;++i)for(size_t j=0;j<np;++j){
    double b=-c->ionosphere_radius_m*c->ionosphere_radius_m*sin(acos(-1.)*i/(nt-1))*(fac[i*np+j]-stats->removed_mean_fac_a_m2);
    res[reduced_index(np,i,j)]-=b;bn+=b*b;
  }
  for(size_t k=0;k<n;++k)rn+=res[k]*res[k];
  double tolerance=fmax(c->absolute_tolerance,c->relative_tolerance*fmax(sqrt(bn),1.));
  fprintf(stderr,"Tensor true residual %.9e target %.9e reported %.9e\n",sqrt(rn),tolerance,stats->solver.final_residual);
  rc=isfinite(rn) && sqrt(rn)<=tolerance ? 0:-1;
  destroy_global_tensor_stencil(&s);free(mem);return rc;
}
