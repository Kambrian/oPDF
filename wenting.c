#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define VecNorm(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#define VecProd(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])

extern double func_model_(double *r, double *vr, double *vt, double *pars);
extern double func_(double *r, double *vr, double *vt);
double dataprob_model(double r, double vr, double vt, double *pars)
{
  return log(func_model_(&r,&vr,&vt, pars));
//   printf("%g,%g,%g,%g,%g,%g: %g,%g,%g; %g\n",pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],r,vr,vt,lnL);
}
double dataprob6d_model(double xv[6], double *pars)
{//dP/d^3xd^3v
  //return log(Prob)
  double *x=xv, *v=xv+3;
  double r=sqrt(VecNorm(x));
  double vr=VecProd(x,v)/r;
  double vt=sqrt(VecNorm(v)-vr*vr);
//   printf("%g\n", log(0.));
  return log(func_model_(&r,&vr,&vt, pars));
}
double dataprob(double a, double b, double c)
{//dP/d^3xd^3v, not dP/drdvrdvt=dP/d^3xd^3v*4*pi*r^2*2*pi*vt
  //return log(Prob)
  return log(func_(&a, &b, &c));
}
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