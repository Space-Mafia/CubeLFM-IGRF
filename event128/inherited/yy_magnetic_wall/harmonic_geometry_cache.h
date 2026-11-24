#ifndef GAMERA_HARMONIC_GEOMETRY_CACHE_H
#define GAMERA_HARMONIC_GEOMETRY_CACHE_H

#include <pthread.h>
#include <stddef.h>
#include <stdlib.h>

/* Two independent fixed-size, zero-initialized tables per live worker. The fast
 * TLS pointer avoids pthread_getspecific in the field-evaluation hot loop.
 * The pthread key owns the allocation: libgomp may retire workers when the
 * root switches between MHD and MI/MKL team sizes. A raw TLS heap pointer
 * alone would lose, not free, that allocation when a worker exits.
 * Slot 0 owns the geometry table, slot 1 the exact point/time vector cache.
 * No cache contents, keys, capacity, or numerical operations change. */
static pthread_key_t gamera_geometry_key[2];
static pthread_once_t gamera_geometry_once = PTHREAD_ONCE_INIT;
static int gamera_geometry_key_ready;

static void gamera_geometry_key_init(void)
{
    if (pthread_key_create(&gamera_geometry_key[0], free) != 0) return;
    if (pthread_key_create(&gamera_geometry_key[1], free) != 0) {
        pthread_key_delete(gamera_geometry_key[0]);
        return;
    }
    gamera_geometry_key_ready = 1;
}

static inline void *gamera_harmonic_thread_cache(unsigned slot, size_t bytes)
{
    static __thread void *table[2];
    static __thread size_t table_bytes[2];
    if (slot >= 2) return NULL;
    if (table[slot]) return bytes == table_bytes[slot] ? table[slot] : NULL;
    if (pthread_once(&gamera_geometry_once, gamera_geometry_key_init) != 0 ||
        !gamera_geometry_key_ready) return NULL;
    void *candidate = calloc(1, bytes);
    if (!candidate) return NULL;
    if (pthread_setspecific(gamera_geometry_key[slot], candidate) != 0) {
        free(candidate);
        return NULL;
    }
    table_bytes[slot] = bytes;
    table[slot] = candidate;
    return table[slot];
}
#endif
