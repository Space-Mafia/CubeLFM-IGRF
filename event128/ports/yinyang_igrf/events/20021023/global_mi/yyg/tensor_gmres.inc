/* Bounded restarted right-preconditioned GMRES fallback for nonsymmetric
 * Hall systems. A and b are unchanged; acceptance uses ||b-Ax||, not the
 * preconditioned norm. The constant-potential nullspace is projected out. */
static int tensor_gmres(gamera_mi_global_tensor_stencil *a,
                        gamera_mi_global_tensor_stencil *p,
                        const double *b, double *x, size_t n,
                        double target, int budget, double *norm) {
  const int m=60;
  double *mem=calloc((2*m+4)*n + (m+1)*m + 4*m+1,sizeof(double));
  if(!mem)return -1;
  double *v=mem,*z=v+(m+1)*n,*r=z+m*n,*w=r+n;
  double *h=w+n,*cs=h+(m+1)*m,*sn=cs+m,*g=sn+m,*y=g+m+1;
  int used=0;
  while(used<budget) {
    apply_global_tensor_stencil(a,x,r);
    for(size_t k=0;k<n;++k)r[k]=b[k]-r[k];
    *norm=sqrt(dot(r,r,n));
    if(!isfinite(*norm) || *norm<=target)break;
    project_zero_sum(r,n);
    double beta=sqrt(dot(r,r,n));
    if(!(beta>DBL_MIN) || !isfinite(beta))break;
    for(size_t k=0;k<n;++k)v[k]=r[k]/beta;
    memset(h,0,(m+1)*m*sizeof(double));
    memset(g,0,(m+1)*sizeof(double));g[0]=beta;
    int q=0;
    for(int j=0;j<m && used<budget;++j) {
      apply_global_tensor_preconditioner(p,v+j*n,z+j*n);
      project_zero_sum(z+j*n,n);
      apply_global_tensor_stencil(a,z+j*n,w);
      /* Twice modified Gram-Schmidt avoids false convergence for the very
       * anisotropic equatorial operator. */
      for(int pass=0;pass<2;++pass)for(int i=0;i<=j;++i) {
        double t=dot(v+i*n,w,n);h[i+(m+1)*j]+=t;
        for(size_t k=0;k<n;++k)w[k]-=t*v[i*n+k];
      }
      double hn=sqrt(dot(w,w,n));h[j+1+(m+1)*j]=hn;
      if(hn>DBL_MIN && isfinite(hn))
        for(size_t k=0;k<n;++k)v[(j+1)*n+k]=w[k]/hn;
      for(int i=0;i<j;++i) {
        double t=cs[i]*h[i+(m+1)*j]+sn[i]*h[i+1+(m+1)*j];
        h[i+1+(m+1)*j]=-sn[i]*h[i+(m+1)*j]+cs[i]*h[i+1+(m+1)*j];
        h[i+(m+1)*j]=t;
      }
      double den=hypot(h[j+(m+1)*j],h[j+1+(m+1)*j]);
      if(!(den>DBL_MIN) || !isfinite(den))break;
      cs[j]=h[j+(m+1)*j]/den;sn[j]=h[j+1+(m+1)*j]/den;
      h[j+(m+1)*j]=den;h[j+1+(m+1)*j]=0.;
      g[j+1]=-sn[j]*g[j];g[j]*=cs[j];
      ++used;q=j+1;
      if(fabs(g[j+1])<=target || hn<=DBL_MIN)break;
    }
    if(q==0)break;
    for(int i=q-1;i>=0;--i) {
      double t=g[i];for(int j=i+1;j<q;++j)t-=h[i+(m+1)*j]*y[j];
      y[i]=t/h[i+(m+1)*i];
    }
    for(int i=0;i<q;++i)for(size_t k=0;k<n;++k)x[k]+=y[i]*z[i*n+k];
    project_zero_sum(x,n);
  }
  apply_global_tensor_stencil(a,x,r);
  for(size_t k=0;k<n;++k)r[k]=b[k]-r[k];
  *norm=sqrt(dot(r,r,n));free(mem);
  return used;
}
