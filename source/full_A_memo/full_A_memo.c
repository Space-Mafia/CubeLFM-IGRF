/* Same qualified bounded bit-key memo policy, adapted value3/raw CGS callback.
 * Setter always invalidates the generation, including same-time reconfiguration.
 * Original field evaluator is linked unchanged; no rounding/interpolation. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#define CACHE_SIZE 16384u
struct entry { uint64_t key[4],generation; double value[3]; };
static struct entry *cache;
static uint64_t generation=1,hits,misses,calls,resets;
static int allocation_attempted,registered;
extern void cube_yy_hydro_potential_uncached(const double*,const double*,double*);
static void stats(void){const char*p=getenv("CUBE_HYDRO_A_MEMO_STATS");if(!p||!*p)return;FILE*f=fopen(p,"wx");if(!f)return;fprintf(f,"{\"calls\":%llu,\"hits\":%llu,\"misses\":%llu,\"resets\":%llu,\"allocated\":%d}\n",(unsigned long long)calls,(unsigned long long)hits,(unsigned long long)misses,(unsigned long long)resets,cache!=NULL);fclose(f);}
void cube_yy_hydro_potential_reset(void){++resets;if(++generation==0){if(cache)memset(cache,0,CACHE_SIZE*sizeof(*cache));generation=1;}}
void cube_yy_hydro_potential_memo(const double*x,double t,double*b){
 ++calls;if(!registered){atexit(stats);registered=1;}
 if(!isfinite(x[0])||!isfinite(x[1])||!isfinite(x[2])||!isfinite(t)){cube_yy_hydro_potential_uncached(x,&t,b);return;}
 if(!allocation_attempted){cache=calloc(CACHE_SIZE,sizeof(*cache));allocation_attempted=1;}
 if(!cache){++misses;cube_yy_hydro_potential_uncached(x,&t,b);return;}
 uint64_t key[4];memcpy(key,x,24);memcpy(key+3,&t,8);uint64_t h=UINT64_C(1469598103934665603);
 for(int i=0;i<4;++i){h=(h^key[i])*UINT64_C(1099511628211);h^=h>>32;}
 struct entry*e=&cache[h&(CACHE_SIZE-1u)];
 if(e->generation==generation&&memcmp(e->key,key,sizeof(key))==0){++hits;memcpy(b,e->value,24);return;}
 ++misses;cube_yy_hydro_potential_uncached(x,&t,b);
 if(isfinite(b[0])&&isfinite(b[1])&&isfinite(b[2])){memcpy(e->key,key,sizeof(key));memcpy(e->value,b,24);e->generation=generation;}
}
