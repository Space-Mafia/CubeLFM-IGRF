#ifndef CUBE_YY_CT_SOURCE_H
#define CUBE_YY_CT_SOURCE_H
#include "nonorthogonal_grid.h"
typedef struct {
 int n[3];
 double length_m;
 gamera_no_grid grid;
 double *flux[3];
 gamera_no_vec3 *magnetic;
} cube_yy_ct_source;
int cube_yy_ct_source_prepare(int,int,int,double,const double*,void**);
int cube_yy_ct_source_apply(void*,const double*,const double*,double*);
void cube_yy_ct_source_destroy(void*);
#endif
