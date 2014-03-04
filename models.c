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

#define TOL_BIN 1e-2 //bin size relative error
#define TOL_REL 1e-5 //relative tolerance for PERIOD integral
#define MAX_INTVAL 5000

double halo_pot(double r)
{
  double x=r/Halo.Rs;
  return Halo.Pots*log(1+x)/x; //the halo boundary is irrelevant in this problem, since only the potential difference affects the orbit
}
static int is_forbidden(double r, int pid)
{
  if(P[pid].E-P[pid].L2/2./r/r-halo_pot(r)<0) 
    return 1;

  return 0;
}

void solve_radial_limits(int pid)
{//estimate pericenter and apocenter (outer bounds, to serve as integration limits)
  P[pid].rlim[0]=R_MIN;
  P[pid].rlim[1]=R_MAX;
  if(!is_forbidden(P[pid].rlim[0],pid)&!is_forbidden(P[pid].rlim[1],pid)) return; //no need to zoom
  
  double xmid;
  xmid=(P[pid].rlim[0]+P[pid].rlim[1])/2;
  while(is_forbidden(xmid, pid)) //shrink the bin till midpoint is allowed
  {
    if(xmid<P[pid].r) 
      P[pid].rlim[0]=xmid;
    else
      P[pid].rlim[1]=xmid;
    xmid=(P[pid].rlim[0]+P[pid].rlim[1])/2;
  }

  //now good mid
    double lbin[2], rbin[2], dx;
    lbin[0]=P[pid].rlim[0];lbin[1]=xmid;
    rbin[0]=xmid;rbin[1]=P[pid].rlim[1];
    dx=(P[pid].rlim[1]-P[pid].rlim[0])/2;
    while(1)
    {
      dx/=2.;
      
      xmid=(lbin[0]+lbin[1])/2.;
      if(is_forbidden(xmid, pid)) 
	lbin[0]=xmid;
      else
	lbin[1]=xmid;
      
      xmid=(rbin[0]+rbin[1])/2.;
      if(is_forbidden(xmid, pid))
	rbin[1]=xmid;
      else
	rbin[0]=xmid;
      
      if(dx/(rbin[1]-lbin[0])<TOL_BIN) break; 
    }
    P[pid].rlim[0]=lbin[0];
    P[pid].rlim[1]=rbin[1];
}

double vr_inv_part(double r, int pid)
{
  double vr2=2*(P[pid].E-P[pid].L2/2./r/r-halo_pot(r));
  if(vr2<0) return 0.;
  return 1./sqrt(vr2);
}
static double vr_inv_rfunc(double r, void *params)
{ /*--- 1/vel_r ---*/
  return vr_inv_part(r, *(int *)params);
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
void solve_radial_orbit(int pid)
{//find peri(apo)-centers and integrate the period 
  P[pid].E=P[pid].K+halo_pot(P[pid].r);
  solve_radial_limits(pid);
  
  gsl_function F;
  F.function = &vr_inv_rfunc;
  F.params = &pid;
  
  double error;
  //   gsl_integration_qags (&F, xlim[0],xlim[1], 0, TOL_REL, MAX_INTVAL, //3, 
  // 		       GSL_workspace, &result, &error);
  size_t neval;
  gsl_integration_cquad (&F, P[pid].rlim[0],P[pid].rlim[1], 0, TOL_REL, //3, 
			 GSL_workspaceC, &(P[pid].T), &error, &neval);
  //   result=smpintD(&F,xlim[0],xlim[1],TOL_REL); //too slow
  if(P[pid].T<=0) printf("Part %d: r=%g,K=%g, E=%g, L2=%g; T=%g (vt/v=%g)\n",pid, P[pid].r, P[pid].K, P[pid].E, P[pid].L2, P[pid].T, sqrt(P[pid].L2/P[pid].r/P[pid].r/2./P[pid].K));
}

double likelihood(double pars[])
{
  int i,j;
  double lnL=0.,p;
  double z=0., M=pars[0], c=pars[1];
  decode_NFWprof(z,M,c,VIR_C200,&Halo);  
  
  #pragma omp parallel 
  {
    #pragma omp for reduction(+:lnL)
    for(i=0;i<nP;i++)
    {
      if(P[i].r<R_MIN||P[i].r>R_MAX) continue;
      solve_radial_orbit(i);
    }
    #pragma omp for private(i,p,j) reduction(+:lnL)
    for(i=0;i<nP;i++)
    {
      for(p=0,j=0;j<nP;j++)
      {
	if(P[i].r<P[j].rlim[0]||P[i].r>P[j].rlim[1]) continue;
	p+=vr_inv_part(P[i].r,j)/P[j].T;//contribution from j to i;
      }
      if(p<=0)
	printf("id=%d, p=%g\n",i, p);
      lnL+=log(p);
    }
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

 double lnL=likelihood(pars);
  printf("M=%g,c=%g,lnL=%g\n",MPAR,CPAR,lnL);

  free_integration_space();
  
  return 0;
}