/* Smooth MLT-resolved boundary. Harmonics1/2 retain offset and ellipticity;
 * no fixed eccentricity or north/south mirroring is imposed. */
static int oval_fit(size_t n,const double *edge,const double *weight,
                    const double *previous,int initialized,double quiet,
                    double minimum,double maximum,double night_limit,
                    double *out,int *limited) {
  double A[5][6]={{0}},mean=0.;size_t valid=0;
  for(size_t j=0;j<n;j++){
    if(!isfinite(weight[j]) || weight[j]<0)return -1;
    if(weight[j]>0 && !isfinite(edge[j]))return -1;
    if(initialized && !isfinite(previous[j]))return -1;
    double a=2.*DPB_PI*j/n,b[5]={1.,cos(a),sin(a),cos(2*a),sin(2*a)};
    double prior=initialized?previous[j]:quiet;
    double w=weight[j],wp=.025;
    if(w>0)valid++;
    double rhs=w>0?w*edge[j]+wp*prior:wp*prior;
    for(int r=0;r<5;r++){
      A[r][5]+=b[r]*rhs/n;
      for(int c=0;c<5;c++)A[r][c]+=(w+wp)*b[r]*b[c]/n;
    }
  }
  /* Mild curvature penalty; evidence, not a hard-coded shape, determines oval. */
  for(int k=1;k<5;k++)A[k][k]+=.002*(k<3?1.:16.);
  if(valid<5){
    for(size_t j=0;j<n;j++)out[j]=initialized?previous[j]:quiet;
    *limited=0;return 0;
  }
  for(int c=0;c<5;c++){
    int p=c;for(int r=c+1;r<5;r++)if(fabs(A[r][c])>fabs(A[p][c]))p=r;
    if(fabs(A[p][c])<1e-14)return -1;
    for(int k=c;k<6;k++){double t=A[p][k];A[p][k]=A[c][k];A[c][k]=t;}
    double v=A[c][c];for(int k=c;k<6;k++)A[c][k]/=v;
    for(int r=0;r<5;r++)if(r!=c){double f=A[r][c];for(int k=c;k<6;k++)A[r][k]-=f*A[c][k];}
  }
  mean=clamp_value(A[0][5],minimum+1e-6,night_limit-1e-6);
  double amplitude=1.;*limited=0;
  for(size_t j=0;j<n;j++){
    double a=2.*DPB_PI*j/n;
    out[j]=A[1][5]*cos(a)+A[2][5]*sin(a)+A[3][5]*cos(2*a)+A[4][5]*sin(2*a);
    double mlt=24.*j/n,top=(mlt<=4.||mlt>=20.)?night_limit:maximum;
    if(out[j]>0)amplitude=fmin(amplitude,(top-mean)/out[j]);
    if(out[j]<0)amplitude=fmin(amplitude,(minimum-mean)/out[j]);
  }
  amplitude=fmax(0.,amplitude);*limited=amplitude<1. || mean!=A[0][5];
  for(size_t j=0;j<n;j++)out[j]=mean+amplitude*out[j];
  return 0;
}
