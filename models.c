#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <sys/times.h>

#include "mymath.h"
#include "cosmology.h"
#include "io.h"
#include "models.h"

struct NFWParZ Halo;
#define R_MIN 30
#define R_MAX 100

#define TOL_BIN 1e-2 //bin size relative error
#define TOL_REL 1e-5 //relative tolerance for PERIOD integral
#define MAX_INTVAL 5000

double halo_pot(double r)
{
  double x=r/Halo.Rs;
  return Halo.Pots*log(1+x)/x; //the halo boundary is irrelevant in this problem, since only the potential difference affects the orbit
}
static int is_forbidden(double r, void *params)
{
  double *p=params; //p[0]=E,p[1]=L^2
  if(p[0]-p[1]/2./r/r-halo_pot(r)<0) 
    return 1;

  return 0;
}
void zoom_bin(double rref, double x[2], void *params)
{
  if(!is_forbidden(x[0],params)&!is_forbidden(x[1],params)) return; //no need to zoom
  
  double xmid;
  xmid=(x[0]+x[1])/2;
  if(is_forbidden(xmid, params)) //
  {
    if(xmid<rref) 
      x[0]=xmid;
    else
      x[1]=xmid;
    zoom_bin(rref, x, params);
  }
  else
  {//now good mid
    double lbin[2], rbin[2], dx;
    lbin[0]=x[0];lbin[1]=xmid;
    rbin[0]=xmid;rbin[1]=x[1];
    dx=(x[1]-x[0])/2;
    while(1)
    {
      dx/=2.;
      
      xmid=(lbin[0]+lbin[1])/2.;
      if(is_forbidden(xmid, params)) 
	lbin[0]=xmid;
      else
	lbin[1]=xmid;
      
      xmid=(rbin[0]+rbin[1])/2.;
      if(is_forbidden(xmid, params))
	rbin[1]=xmid;
      else
	rbin[0]=xmid;
      
      if(dx/(rbin[1]-lbin[0])<TOL_BIN) break; 
    }
    x[0]=lbin[0];
    x[1]=rbin[1];
  }
}

double vr_inv(double r, void *params)
{ /*--- 1/vel_r ---*/
  double *p=params; //p[0]=E,p[1]=L^2
  double vr2=2*(p[0]-p[1]/2./r/r-halo_pot(r));
  if(vr2<0) return 0.;
  return 1./sqrt(vr2);
}
// static gsl_integration_workspace * GSL_workspace;
// #pragma omp threadprivate(GSL_workspaceC)
static gsl_integration_cquad_workspace * GSL_workspaceC;
#pragma omp threadprivate(GSL_workspaceC)
void alloc_integration_space()
{//have to allocate two workspaces when evaluating double integral, to avoid entangling inner and outer workspaces.
	#pragma omp parallel
  {//GSL_workspace=gsl_integration_workspace_alloc(MAX_INTVAL);
	GSL_workspaceC=gsl_integration_cquad_workspace_alloc(MAX_INTVAL);}
}
void free_integration_space()
{
	#pragma omp parallel
  {//gsl_integration_workspace_free (GSL_workspace);
    gsl_integration_cquad_workspace_free (GSL_workspaceC);}
}
double Period(double r, double K, double L2, double rmin, double rmax)
{
 //integrate the period 
  double p[2];
  p[0]=K+halo_pot(r); //E=K+pot
  p[1]=L2; //L^2
//   printf("%g,%g,%g,%g\n",r,K,p[0],p[1]);
  gsl_function F;
  F.function = &vr_inv;
  F.params = p;
  double xlim[2];
  xlim[0]=rmin;xlim[1]=rmax;
  zoom_bin(r,xlim,p);
//   printf("bin:[%g,%g]\n", xlim[0], xlim[1]);
//   double x;
//   for(x=xlim[0];x<xlim[1];x+=(xlim[1]-xlim[0])/20)
//     printf("r=%g,1/v=%g\n",x,vr_inv(x,p));
  double result,error;
//   gsl_integration_qags (&F, xlim[0],xlim[1], 0, TOL_REL, MAX_INTVAL, //3, 
// 		       GSL_workspace, &result, &error);
  size_t neval;
  gsl_integration_cquad (&F, xlim[0],xlim[1], 0, TOL_REL, //3, 
		       GSL_workspaceC, &result, &error, &neval);
  //   result=smpintD(&F,xlim[0],xlim[1],TOL_REL); //too slow
  if(result<=0) printf("%g,%g(%g),%g: %g (vt/v=%g)\n", r, K, p[0], L2, result, sqrt(L2/r/r/2./K));
  return result; 
}

double likelihood(double pars[])
{
  int i;
  double lnL=0.;
  double z=0., M=pars[0], c=pars[1];
  decode_NFWprof(z,M,c,VIR_C200,&Halo);  
  
  #pragma omp parallel for reduction(+:lnL)
  for(i=0;i<nP;i++)
  {
    if(P[i].r<R_MIN||P[i].r>R_MAX) continue;
    double T=Period(P[i].r,P[i].K,P[i].L2, R_MIN,R_MAX);
    if(T<=0) 
    {
      printf("%g,%g,%g\n", P[i].r, P[i].K+halo_pot(P[i].r), P[i].L2);
      printf("%d:T=%g,[%g,%g,%g],[%g,%g,%g]\n",i,T,P[i].x[0],P[i].x[1],P[i].x[2],P[i].v[0],P[i].v[1],P[i].v[2]); 
      printf("rs=%g, pots=%g\n",Halo.Rs,Halo.Pots);
      exit(1);
    };
    lnL+=-log(T);
//     lnL+=-log(Period(P[i].r,P[i].K,P[i].L2)); //scale the period with L/K?
  }
  
  return lnL;
}

#define MPAR pars[0]
#define CPAR pars[1]
int main(int argc, char **argv)
{
  char datafile[1024]=ROOTDIR"/data/mockhalo.hdf5";
  double pars[NUM_PAR_MAX]={2e2,15};
  if(argc==3)
  {
    MPAR=atof(argv[1]);
    CPAR=atof(argv[2]);
  }
  
  load_data(datafile);
  alloc_integration_space();
//   decode_NFWprof(0,MPAR,CPAR,VIR_C200,&Halo);
//   int i=618;
//   printf("T=%g\n",Period(P[i].r,P[i].K,P[i].L2,63,63.3));
 double lnL=likelihood(pars);
  printf("M=%g,c=%g,lnL=%g\n",MPAR,CPAR,lnL);

  free_integration_space();
  
  return 0;
}