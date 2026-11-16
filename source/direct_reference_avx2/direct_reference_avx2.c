/* Same ordered21 FP64 products, fixed independent32 output lanes.
 * Write directly in the existing Fortran destination's contiguous layout.
 * No horizontal reduction, interpolation, phase or coefficient change. */
#include <immintrin.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
int cube_yy_direct_reference_evaluate(int count,const double *restrict bank,
    const int *restrict active,const double *restrict weight,double *restrict output){
 if(count<0||!bank||!active||!weight||!output)return 1;
 for(int k=0;k<21;++k)if(!isfinite(weight[k]))return 1;
 const int blocks=(count+31)/32;
 for(int b=0;b<blocks;++b){
  const int remaining=count-b*32,n=remaining<32?remaining:32;
  if(!active[b]){memset(output+(size_t)b*32,0,(size_t)n*sizeof(double));continue;}
  __m256d s0=_mm256_setzero_pd(),s1=s0,s2=s0,s3=s0,s4=s0,s5=s0,s6=s0,s7=s0;
  const double *a=bank+(size_t)b*32*21;
  for(int k=0;k<21;++k){
   const __m256d w=_mm256_set1_pd(weight[k]);const double *p=a+k*32;
   s0=_mm256_add_pd(s0,_mm256_mul_pd(_mm256_loadu_pd(p),w));
   s1=_mm256_add_pd(s1,_mm256_mul_pd(_mm256_loadu_pd(p+4),w));
   s2=_mm256_add_pd(s2,_mm256_mul_pd(_mm256_loadu_pd(p+8),w));
   s3=_mm256_add_pd(s3,_mm256_mul_pd(_mm256_loadu_pd(p+12),w));
   s4=_mm256_add_pd(s4,_mm256_mul_pd(_mm256_loadu_pd(p+16),w));
   s5=_mm256_add_pd(s5,_mm256_mul_pd(_mm256_loadu_pd(p+20),w));
   s6=_mm256_add_pd(s6,_mm256_mul_pd(_mm256_loadu_pd(p+24),w));
   s7=_mm256_add_pd(s7,_mm256_mul_pd(_mm256_loadu_pd(p+28),w));
  }
  double tail[32],*r=n==32?output+(size_t)b*32:tail;
  _mm256_storeu_pd(r,s0);_mm256_storeu_pd(r+4,s1);_mm256_storeu_pd(r+8,s2);_mm256_storeu_pd(r+12,s3);
  _mm256_storeu_pd(r+16,s4);_mm256_storeu_pd(r+20,s5);_mm256_storeu_pd(r+24,s6);_mm256_storeu_pd(r+28,s7);
  if(n!=32)memcpy(output+(size_t)b*32,tail,(size_t)n*sizeof(double));
 }
 return 0;
}
