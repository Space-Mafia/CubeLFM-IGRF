/* SOURCE-ONLY proposal. Integer classification, no floating arithmetic. */
#include <float.h>
#include <immintrin.h>
#include <stdint.h>
#include <string.h>

_Static_assert(sizeof(double)==8 && FLT_RADIX==2 && DBL_MANT_DIG==53 && DBL_MAX_EXP==1024,
               "The original REAL64 input contract requires IEEE binary64");
_Static_assert(sizeof(int64_t)==8, "Fortran SIZE(kind=8) count must be int64");
#ifndef __STDC_IEC_559__
#error "This fast predicate requires the original IEEE floating representation"
#endif

/* Existing Fortran ABI: INTEGER function; address of array; int64 count byref.
 * A zero return does not replace the original invalid-input behavior: the
 * wrapper executes its unchanged IEEE rejection expression on that path.
 */
int lfm_cube_all_finite_fast_(const double *array,const int64_t *count) {
    const uint64_t exponent=UINT64_C(0x7ff0000000000000);
    if(*count<0) return 0;
    const uint64_t n=(uint64_t)*count;
    const unsigned char *p=(const unsigned char *)array;
    const __m256i mask=_mm256_set1_epi64x((long long)exponent);
    uint64_t i=0;
    while(n-i>=4) {
        __m256i bits;
        memcpy(&bits,p+i*8,sizeof(bits));
        const __m256i rejected=_mm256_cmpeq_epi64(_mm256_and_si256(bits,mask),mask);
        if(!_mm256_testz_si256(rejected,rejected)) return 0;
        i+=4;
    }
    for(;i<n;i++) {
        uint64_t bits;
        memcpy(&bits,p+i*8,sizeof(bits));
        if((bits&exponent)==exponent) return 0;
    }
    return 1;
}
