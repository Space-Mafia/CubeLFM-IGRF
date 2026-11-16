/* Memoize only a previously PASSED unchanged geometric predicate.
 * Every lookup compares every FRAME byte; mutation, signed-zero changes,
 * new shape, eviction or allocation failure use the original Fortran checks.
 * This never caches a fluid, magnetic, positivity or other dynamic check.
 */
#include <fenv.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
_Static_assert(sizeof(double)==8&&sizeof(int64_t)==8,"Original REAL8/SIZEkind8 ABI");
#define FRAME_SLOTS 256u
struct slot { const void *address; int64_t count; unsigned char *words; int rounding; };
static struct slot slots[FRAME_SLOTS];
static int64_t hits,misses,records;
static int valid(const void *p,const int64_t *n){return p&&n&&*n>0&&(uint64_t)*n<=SIZE_MAX/sizeof(double);}
static unsigned index_of(const void *p){uintptr_t x=(uintptr_t)p;x^=x>>17;x^=x>>9;return(unsigned)x&(FRAME_SLOTS-1u);}
static struct slot *locate(const void *p){
 unsigned h=index_of(p);
 for(unsigned i=0;i<FRAME_SLOTS;++i){struct slot *s=&slots[(h+i)&(FRAME_SLOTS-1u)];if(!s->address||s->address==p)return s;}
 return &slots[h];
}
int lfm_cube_frame_validation_hit_(const double *frame,const int64_t *count){
 if(!valid(frame,count)){++misses;return 0;}
 struct slot *s=locate(frame);
 if(s->address==frame&&s->count==*count&&s->rounding==fegetround()&&s->words&&
    memcmp(s->words,frame,(size_t)*count*sizeof(double))==0){++hits;return 1;}
 ++misses;return 0;
}
void lfm_cube_frame_validation_record_(const double *frame,const int64_t *count){
 if(!valid(frame,count))return;
 struct slot *s=locate(frame);const size_t bytes=(size_t)*count*sizeof(double);
 if(!s->words||s->count!=*count){
  unsigned char *copy=malloc(bytes);
  if(!copy)return; /* The just-completed original checks remain authoritative. */
  free(s->words);s->words=copy;
 }
 memcpy(s->words,frame,bytes);s->address=frame;s->count=*count;s->rounding=fegetround();++records;
}
void lfm_cube_frame_validation_statistics_(int64_t counts[3]){counts[0]=hits;counts[1]=misses;counts[2]=records;}
