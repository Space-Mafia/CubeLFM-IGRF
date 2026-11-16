#ifndef CUBE_YY_MAGNETIC_WALL_H
#define CUBE_YY_MAGNETIC_WALL_H
#include <stddef.h>
#include "nonorthogonal_geometry.h"
typedef struct {
  gamera_no_vec3 point;
  double vector[3], norm, scalar;
} cube_yy_wall_face;
typedef struct {
  int ng,nj,nk;
  double length_cm;
  cube_yy_wall_face *face;
} cube_yy_magnetic_wall;
int cube_yy_magnetic_wall_create(int,int,int,const double*,const double*,const double*,void**);
int cube_yy_magnetic_wall_apply(void*,int,double,double,double,double*);
void cube_yy_magnetic_wall_destroy(void*);
int cube_yy_magnetic_wall_field(const double[3],double,double[3]);
/* Fortran order: cached radial 1-ng..1, angular 1-ng..n+ng, direction. */
static inline size_t cube_yy_wall_index(const cube_yy_magnetic_wall *w,int i,int j,int k,int d) {
  return i+(w->ng+1)*(j+(w->nj+2*w->ng)*(k+(w->nk+2*w->ng)*d));
}
#endif
