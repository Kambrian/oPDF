#include <math.h>
#include <stdlib.h>

extern double func_(double *r, double *vr, double *vt);
double dataprob(double a, double b, double c)
{//dP/d^3xd^3v, not dP/drdvrdvt=dP/d^3xd^3v*4*pi*r^2*2*pi*vt
  //return log(Prob)
  return log(func_(&a, &b, &c));
}
#define VecNorm(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#define VecProd(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
double dataprob6d(double xv[6])
{//dP/d^3xd^3v
  //return log(Prob)
  double *x=xv, *v=xv+3;
  double r=sqrt(VecNorm(x));
  double vr=VecProd(x,v)/r;
  double vt=sqrt(VecNorm(v)-vr*vr);
//   printf("%g\n", log(0.));
  return log(func_(&r,&vr,&vt));
}