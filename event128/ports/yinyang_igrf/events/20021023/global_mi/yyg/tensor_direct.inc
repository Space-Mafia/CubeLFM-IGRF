/* Optional MKL sparse-LU backend. No stencil, tensor, source, tolerance or
 * physics changes. The independent checked_tensor.c acceptance is retained.
 * Each changing conductance is refactorized: no stale numerical factors.
 * No shared mutable cache; root-only caller, MKL uses its allotted threads. */
#include <mkl.h>
#include <limits.h>
#include <stdio.h>

static void tensor_area_gauge(const gamera_mi_global_solver_config *c,double *x){
  size_t np=c->longitude_count,nt=c->colatitude_count;
  double sum=0.,weight=0.,dtheta=acos(-1.)/(nt-1);
  for(size_t i=1;i+1<nt;++i)for(size_t j=0;j<np;++j){
    double w=sin(i*dtheta);sum+=w*x[reduced_index(np,i,j)];weight+=w;
  }
  double mean=sum/weight;
  for(size_t k=0;k<2+np*(nt-2);++k)x[k]-=mean;
}

static int tensor_sparse_direct(gamera_mi_global_tensor_stencil *a,
    const double *b,double *x,size_t count,double target,int budget,
    int *steps,double *norm){
  size_t np=a->config->longitude_count,nt=a->config->colatitude_count;
  size_t n=count-1,capacity=5*n+np;
  *steps=0;*norm=INFINITY;
  if(n>INT_MAX || capacity>INT_MAX || !isfinite(target))return -1;
  MKL_INT *ia=calloc(n+1,sizeof(*ia)),*ja=calloc(capacity,sizeof(*ja));
  double *val=calloc(capacity,sizeof(*val));
  double *mem=calloc(3*count,sizeof(double));
  if(!ia || !ja || !val || !mem){free(ia);free(ja);free(val);free(mem);return -1;}
  double *rhs=mem,*delta=mem+count,*res=mem+2*count;
  size_t z=0;
  for(size_t k=1;k<count;++k){
    ia[k-1]=(MKL_INT)z;
    if(k==count-1){
      for(size_t j=0;j<np;++j){
        ja[z]=(MKL_INT)(reduced_index(np,nt-2,j)-1);val[z++]=-a->pole_edge_south[j];
      }
      ja[z]=(MKL_INT)(k-1);val[z++]=a->diagonal[k];
    }else{
      size_t i=(k-1)/np+1,j=(k-1)%np;
      size_t cols[5]={k,i==1?0:k-np,i==nt-2?count-1:k+np,
        j==0?k+np-1:k-1,j==np-1?k-np+1:k+1};
      double vv[5]={a->diagonal[k],a->lower[k],a->upper[k],
        a->longitude_lower[k],a->longitude_upper[k]};
      for(int q=1;q<5;++q)for(int r=q;r>0 && cols[r]<cols[r-1];--r){
        size_t c=cols[r];cols[r]=cols[r-1];cols[r-1]=c;
        double v=vv[r];vv[r]=vv[r-1];vv[r-1]=v;
      }
      for(int q=0;q<5;++q)if(cols[q]!=0){ja[z]=(MKL_INT)(cols[q]-1);val[z++]=vv[q];}
    }
  }
  ia[n]=(MKL_INT)z;
  void *pt[64]={0};MKL_INT iparm[64]={0},mtype=11;
  pardisoinit(pt,&mtype,iparm);
  iparm[34]=1; /* zero-based CSR */
  iparm[26]=1; /* matrix checker */
  iparm[7]=4;  /* internal iterative refinement, not acceptance criterion */
  MKL_INT maxfct=1,mnum=1,phase=12,nn=(MKL_INT)n,nrhs=1,msg=0,error=0,dummy=0;
  pardiso(pt,&maxfct,&mnum,&mtype,&phase,&nn,val,ia,ja,&dummy,&nrhs,iparm,&msg,rhs,delta,&error);
  int rc=-1;
  if(error){fprintf(stderr,"MI sparse factor error %d\n",(int)error);goto cleanup;}
  memset(x,0,count*sizeof(double));
  /* Explicit residual refinement uses the original full stencil, never a
   * scaled/factorization norm. Small bounded cap avoids endless stagnation. */
  for(int it=0;it<budget && it<16;++it){
    apply_global_tensor_stencil(a,x,res);
    for(size_t k=0;k<count;++k)res[k]=b[k]-res[k];
    *norm=sqrt(dot(res,res,count));
    if(isfinite(*norm) && *norm<=target){rc=0;break;}
    if(!isfinite(*norm))break;
    /* Project only the correction RHS onto the roundoff-compatible range;
     * acceptance above/below always uses the unprojected full residual. */
    project_zero_sum(res,count);
    for(size_t k=1;k<count;++k)rhs[k-1]=res[k];
    phase=33;
    pardiso(pt,&maxfct,&mnum,&mtype,&phase,&nn,val,ia,ja,&dummy,&nrhs,iparm,&msg,rhs,delta,&error);
    if(error){fprintf(stderr,"MI sparse solve error %d\n",(int)error);break;}
    for(size_t k=1;k<count;++k)x[k]+=delta[k-1];
    tensor_area_gauge(a->config,x);
    ++*steps;
    /* Keep north reference zero during correction; output's area-weighted
     * gauge and resulting residual are checked independently by the caller. */
  }
  apply_global_tensor_stencil(a,x,res);
  for(size_t k=0;k<count;++k)res[k]=b[k]-res[k];
  *norm=sqrt(dot(res,res,count));
  rc=isfinite(*norm) && *norm<=target?0:-1;
cleanup:
  phase=-1;error=0;
  pardiso(pt,&maxfct,&mnum,&mtype,&phase,&nn,val,ia,ja,&dummy,&nrhs,iparm,&msg,rhs,delta,&error);
  if(error)rc=-1;
  free(ia);free(ja);free(val);free(mem);return rc;
}
