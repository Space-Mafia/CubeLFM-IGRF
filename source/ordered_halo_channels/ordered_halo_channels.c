/* Original donor order; independent primitive channels in AVX2 lanes.
 * Bounded modes and unsupported/invalid plan metadata use the unchanged
 * Fortran path. Original finite and positive donor screens run before this.
 * No reconstruction, clipping, floors, interpolation weights or MPI changes. */
#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
_Static_assert(sizeof(double)==8&&sizeof(int)==4,"Cube REAL64/INTEGER C ABI");
int cube_yy_halo_ordered8(int nt,int nd,int ldr,int ldw,const int *terms,
 const int *reference,const double *weight,const double *samples,const int *positive,
 double *result,int *limited){
 if(nt<0||nd<0||ldr<1||ldw<1||!terms||!reference||!weight||!samples||!positive||!result||!limited)return 2;
 *limited=0;
 /* Validate every index/weight before any gather. Unsupported metadata
  * returns to the original routine, preserving its original failure policy. */
 for(int h=0;h<nt;++h){
  int d=terms[h];if(d<1||d>ldr||d>ldw)return 2;
  for(int k=0;k<d;++k)if(reference[(size_t)h*ldr+k]<1||reference[(size_t)h*ldr+k]>nd||!isfinite(weight[(size_t)h*ldw+k]))return 2;
 }
 for(int h=0;h<nt;++h){
  const int d=terms[h];const int *r=reference+(size_t)h*ldr;
  const double *w=weight+(size_t)h*ldw;double *out=result+(size_t)h*8;
  if(d==1){memcpy(out,samples+(size_t)(r[0]-1)*8,8*sizeof(double));continue;}
  __m256d lo=_mm256_setzero_pd(),hi=lo;
  for(int k=0;k<d;++k){
   const double *p=samples+(size_t)(r[k]-1)*8;const __m256d a=_mm256_set1_pd(w[k]);
   lo=_mm256_add_pd(lo,_mm256_mul_pd(a,_mm256_loadu_pd(p)));
   hi=_mm256_add_pd(hi,_mm256_mul_pd(a,_mm256_loadu_pd(p+4)));
  }
  _mm256_storeu_pd(out,lo);_mm256_storeu_pd(out+4,hi);
  for(int z=0;z<8;++z)if(positive[z]&&out[z]<=0.){
   double minimum=samples[(size_t)(r[0]-1)*8+z];
   for(int k=1;k<d;++k){double x=samples[(size_t)(r[k]-1)*8+z];if(x<minimum)minimum=x;}
   out[z]=minimum;++*limited;
  }
 }
 return 0;
}
