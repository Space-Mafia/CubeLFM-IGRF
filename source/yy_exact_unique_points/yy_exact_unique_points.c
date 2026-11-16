/* Exact coordinate grouping only: same bit key as the qualified point map.
 * Hash collisions compare all24 bytes. Cyclic work ordering balances ranks.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
static uint64_t rotate(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static uint64_t hash(const uint64_t key[3]){
    uint64_t x=key[0]^rotate(key[1],21)^rotate(key[2],43);return x^(x>>29);
}
int cube_yy_exact_unique_points(const int64_t *keys,int n,int ranks,int *map,int64_t *unique){
    if(n<0 || ranks<1 || (n && (!keys||!map||!unique)))return -1;
    if(!n)return 0;
    size_t capacity=1;
    while(capacity<2*(size_t)n){if(capacity>SIZE_MAX/2)return -1;capacity*=2;}
    int *slots=calloc(capacity,sizeof(int));int *order=malloc((size_t)n*sizeof(int));
    int64_t *work=malloc((size_t)n*24);
    if(!slots||!order||!work){free(slots);free(order);free(work);return -1;}
    int count=0;
    for(int i=0;i<n;++i){
        uint64_t key[3];memcpy(key,keys+3*(size_t)i,24);
        size_t slot=hash(key)&(capacity-1);
        while(slots[slot] && memcmp(unique+3*(size_t)(slots[slot]-1),key,24))slot=(slot+1)&(capacity-1);
        if(!slots[slot]){memcpy(unique+3*(size_t)count,key,24);slots[slot]=++count;}
        map[i]=slots[slot];
    }
    int next=0;
    for(int rank=0;rank<ranks;++rank)
        for(int j=rank;j<count;j+=ranks){
            order[j]=next;memcpy(work+3*(size_t)next,unique+3*(size_t)j,24);++next;
        }
    for(int i=0;i<n;++i)map[i]=order[map[i]-1]+1;
    memcpy(unique,work,(size_t)count*24);
    free(slots);free(order);free(work);return count;
}
