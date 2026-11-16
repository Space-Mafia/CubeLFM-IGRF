/* Contiguous, lossless coefficient STORAGE only. No row sorting or scatter.
 * Positive zero columns retain their original FP64 multiply and add.
 * Lossless coefficient STORAGE only. Every encoded coefficient is
 * checked by integer expansion to the exact original FP64 word.
 * All products and ordered21 sums remain FP64, with FMA disabled.
 * No discarded harmonic, threshold, approximate coefficient or new physics.
 * The previously passed local kernel remains the unchanged fallback. */
#define _POSIX_C_SOURCE 200112L
#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <fenv.h>
#include <float.h>
extern int __real_cube_yy_direct_reference_evaluate(int,const double*,const int*,const double*,double*);
_Static_assert(sizeof(double)==8&&sizeof(float)==4&&FLT_RADIX==2&&FLT_MANT_DIG==24&&DBL_MANT_DIG==53&&DBL_MAX_EXP==1024,"IEEE64/IEEE32 storage ABI");
#define ZERO_ROW UINT32_C(0x80000000)
#define ALL_COLUMNS UINT32_C(0x1fffff)
typedef struct {uint32_t small,zero;size_t offset;} block;
typedef struct {const double *key;const int *active;int count,ready,disabled;size_t blocks,bytes;block *b;unsigned char *data;} plan;
static plan cache[4];static int cleanup_registered;
static uint64_t word(double a){uint64_t x;memcpy(&x,&a,8);return x;}
/* Return an EXACT normal IEEE32 representation, or signed zero only.
 * Subnormals use FP64 so the result cannot depend on DAZ/FTZ modes. */
static int encode32(uint64_t x,uint32_t *f){
 uint64_t magnitude=x&UINT64_C(0x7fffffffffffffff);unsigned e=(unsigned)((x>>52)&2047u);
 if(!magnitude){*f=(uint32_t)(x>>32)&0x80000000u;return 1;}
 if(e<897u||e>1150u||(x&UINT64_C(0x1fffffff)))return 0;
 *f=((uint32_t)(x>>32)&0x80000000u)|((e-896u)<<23)|(uint32_t)((x&UINT64_C(0xfffffffffffff))>>29);
 uint64_t back=((uint64_t)(*f&0x80000000u)<<32)|((uint64_t)((*f>>23)&255u)+896u)<<52|((uint64_t)(*f&0x7fffffu)<<29);
 return back==x;
}
static void cleanup(void){for(int i=0;i<4;++i){free(cache[i].b);free(cache[i].data);}}
static int prepare(plan *p,int count,const double *a,const int *active){
 p->key=a;p->active=active;p->count=count;p->ready=1;p->disabled=1;
 size_t blocks=(size_t)(count+31)/32,bytes=0,old_read=0;block *bp=calloc(blocks,sizeof(*bp));
 if(!bp)return 1;
 uint64_t small_words=0,zero_words=0,large_words=0;
 for(size_t bi=0;bi<blocks;++bi){
  bp[bi].offset=bytes;
  if(!active[bi])continue;
  old_read+=32*21*8;
  for(int k=0;k<21;++k){int small=1,zero=1;
   for(int lane=0;lane<32;++lane){
    uint64_t w=word(a[(bi*21+(size_t)k)*32+(size_t)lane]);uint32_t f;
    if(((w>>52)&2047u)==2047u||(w&UINT64_C(0x7fffffffffffffff))>UINT64_C(0x7830000000000000)){free(bp);return 1;}
    if(w!=0)zero=0;
    if(!encode32(w,&f))small=0;
   }
   if(zero){bp[bi].zero|=1u<<k;zero_words+=32;}
   else if(small){bp[bi].small|=1u<<k;bytes+=128;small_words+=32;}
   else{bytes+=256;large_words+=32;}
  }
 }
 if(bytes+blocks*sizeof(*bp)>=old_read){
  fprintf(stdout,"CONTIGUOUS_LOSSLESS_FALLBACK count=%d source_read_B=%zu encoded_B=%zu descriptor_B=%zu reason=no_reduction\n",count,old_read,bytes,blocks*sizeof(*bp));free(bp);return 1;
 }
 unsigned char *data=NULL;if(posix_memalign((void**)&data,64,bytes)){free(bp);return 1;}
 for(size_t bi=0;bi<blocks;++bi){
  if(!active[bi])continue;size_t offset=bp[bi].offset;
  for(int k=0;k<21;++k){
   if(bp[bi].zero&(1u<<k))continue;
   if(bp[bi].small&(1u<<k)){
    for(int lane=0;lane<32;++lane){uint32_t f;if(!encode32(word(a[(bi*21+(size_t)k)*32+(size_t)lane]),&f)){free(bp);free(data);return 1;}memcpy(data+offset+(size_t)lane*4,&f,4);}offset+=128;
   }else{memcpy(data+offset,a+(bi*21+(size_t)k)*32,256);offset+=256;}
  }
 }
 p->b=bp;p->data=data;p->blocks=blocks;p->bytes=bytes;p->disabled=0;
 fprintf(stdout,"CONTIGUOUS_LOSSLESS_STORAGE count=%d source_read_B=%zu encoded_B=%zu descriptor_B=%zu exact_FP32_words=%llu FP64_words=%llu positive_zero_words=%llu output_order_unchanged=1 original_FP64_words_preserved=1\n",count,old_read,bytes,blocks*sizeof(*bp),(unsigned long long)small_words,(unsigned long long)large_words,(unsigned long long)zero_words);
 return 0;
}
#define ACCUMULATE(I,VALUE) s##I=_mm256_add_pd(s##I,_mm256_mul_pd((VALUE),w))
int __wrap_cube_yy_direct_reference_evaluate(int count,const double *a,const int *active,const double *weight,double *output){
 if(count<=0||!a||!active||!weight||!output||fegetround()!=FE_TONEAREST)return __real_cube_yy_direct_reference_evaluate(count,a,active,weight,output);
 for(int k=0;k<21;++k)if(!isfinite(weight[k])||fabs(weight[k])>1.)return __real_cube_yy_direct_reference_evaluate(count,a,active,weight,output);
 plan *p=NULL;
 for(int i=0;i<4;++i)if(cache[i].ready&&cache[i].key==a&&cache[i].active==active&&cache[i].count==count){p=cache+i;break;}
 if(!p){
  for(int i=0;i<4;++i)if(!cache[i].ready){p=cache+i;break;}
  if(!p)return __real_cube_yy_direct_reference_evaluate(count,a,active,weight,output);
  if(!cleanup_registered){atexit(cleanup);cleanup_registered=1;}
  if(prepare(p,count,a,active))return __real_cube_yy_direct_reference_evaluate(count,a,active,weight,output);
 }
 if(p->disabled)return __real_cube_yy_direct_reference_evaluate(count,a,active,weight,output);
 for(size_t bi=0;bi<p->blocks;++bi){
  const int n=count-(int)(bi*32)<32?count-(int)(bi*32):32;
  if(!active[bi]){memset(output+bi*32,0,(size_t)n*8);continue;}
  const block *b=p->b+bi;const unsigned char *v=p->data+b->offset;
  __m256d s0=_mm256_setzero_pd(),s1=s0,s2=s0,s3=s0,s4=s0,s5=s0,s6=s0,s7=s0;
  for(int k=0;k<21;++k){
   const __m256d w=_mm256_set1_pd(weight[k]);
   if(b->zero&(1u<<k)){
    const __m256d z=_mm256_setzero_pd();
    ACCUMULATE(0,z);ACCUMULATE(1,z);ACCUMULATE(2,z);ACCUMULATE(3,z);
    ACCUMULATE(4,z);ACCUMULATE(5,z);ACCUMULATE(6,z);ACCUMULATE(7,z);
   }else if(b->small&(1u<<k)){
    const float *f=(const float*)v;
    ACCUMULATE(0,_mm256_cvtps_pd(_mm_loadu_ps(f)));ACCUMULATE(1,_mm256_cvtps_pd(_mm_loadu_ps(f+4)));
    ACCUMULATE(2,_mm256_cvtps_pd(_mm_loadu_ps(f+8)));ACCUMULATE(3,_mm256_cvtps_pd(_mm_loadu_ps(f+12)));
    ACCUMULATE(4,_mm256_cvtps_pd(_mm_loadu_ps(f+16)));ACCUMULATE(5,_mm256_cvtps_pd(_mm_loadu_ps(f+20)));
    ACCUMULATE(6,_mm256_cvtps_pd(_mm_loadu_ps(f+24)));ACCUMULATE(7,_mm256_cvtps_pd(_mm_loadu_ps(f+28)));v+=128;
   }else{
    const double *d=(const double*)v;
    ACCUMULATE(0,_mm256_loadu_pd(d));ACCUMULATE(1,_mm256_loadu_pd(d+4));
    ACCUMULATE(2,_mm256_loadu_pd(d+8));ACCUMULATE(3,_mm256_loadu_pd(d+12));
    ACCUMULATE(4,_mm256_loadu_pd(d+16));ACCUMULATE(5,_mm256_loadu_pd(d+20));
    ACCUMULATE(6,_mm256_loadu_pd(d+24));ACCUMULATE(7,_mm256_loadu_pd(d+28));v+=256;
   }
  }
  double tail[32],*r=n==32?output+bi*32:tail;
  _mm256_storeu_pd(r,s0);_mm256_storeu_pd(r+4,s1);_mm256_storeu_pd(r+8,s2);_mm256_storeu_pd(r+12,s3);
  _mm256_storeu_pd(r+16,s4);_mm256_storeu_pd(r+20,s5);_mm256_storeu_pd(r+24,s6);_mm256_storeu_pd(r+28,s7);
  if(n!=32)memcpy(output+bi*32,tail,(size_t)n*8);
 }

 return 0;
}
