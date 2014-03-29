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

double halo_pot(double r)
{
  double x=r/Halo.Rs;
  return Halo.Pots*log(1+x)/x; //the halo boundary is irrelevant in this problem, since only the potential difference affects the orbit
}
static int is_forbidden(double r, int pid)
{//according to current potential and E,L,r
  if(P[pid].E-P[pid].L2/2./r/r-halo_pot(r)<0) 
    return 1;

  return 0;
}

void solve_radial_limits(int pid)
{//estimate pericenter and apocenter (outer bounds, to serve as integration limits), 
  //according to current potential and E,L,r
  P[pid].rlim[0]=R_MIN;
  P[pid].rlim[1]=R_MAX;
  if(is_forbidden(P[pid].r,pid)) //forbidden due to large difference between current pot and initial pot
  {
    P[pid].rlim[0]=P[pid].r;
    P[pid].rlim[1]=P[pid].r;
    return;
  }
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
      
      if(dx/(rbin[1]-lbin[0])<MODEL_TOL_BIN) break; 
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
  {//GSL_workspace=gsl_integration_workspace_alloc(MODEL_MAX_INTVAL);
	GSL_workspaceC=gsl_integration_cquad_workspace_alloc(MODEL_MAX_INTVAL);}
}
void free_integration_space()
{
	#pragma omp parallel
  {//gsl_integration_workspace_free (GSL_workspace);
    gsl_integration_cquad_workspace_free (GSL_workspaceC);}
}
void solve_radial_orbit(int pid)
{//find peri(apo)-centers and integrate the period 
  solve_radial_limits(pid); //according to current potential and E,L,r; E must be initialized with a initial potential
  if(is_forbidden(P[pid].r,pid))
  {
    P[pid].T=1.;  //arbitrary, just to avoid 1/v/T=NaN, since 1/v=0.
    return;
  }
  
  gsl_function F;
  F.function = &vr_inv_rfunc;
  F.params = &pid;
  
  double error;
  //   gsl_integration_qags (&F, xlim[0],xlim[1], 0, MODEL_TOL_REL, MODEL_MAX_INTVAL, //3, 
  // 		       GSL_workspace, &result, &error);
  size_t neval;
  gsl_integration_cquad (&F, P[pid].rlim[0],P[pid].rlim[1], 0, MODEL_TOL_REL, //3, 
			 GSL_workspaceC, &(P[pid].T), &error, &neval);
#if ESTIMATOR==RADIAL_PHASE_ESTIMATOR
  double t;
  gsl_integration_cquad (&F, P[pid].rlim[0],P[pid].r, 0, MODEL_TOL_REL, //3, 
			 GSL_workspaceC, &t, &error, &neval);
  t=t/P[pid].T/2.;
  P[pid].theta=P[pid].vr>=0?t:(1-t);//radial phase
#endif  
  //   result=smpintD(&F,xlim[0],xlim[1],MODEL_TOL_REL); //too slow
//   if(P[pid].T<=0) printf("Part %d (M=%g,c=%g): r=%g,K=%g, E=%g, L2=%g; T=%g (vt/v=%g)\n",pid, Halo.M, Halo.c, P[pid].r, P[pid].K, P[pid].E, P[pid].L2, P[pid].T, sqrt(P[pid].L2/P[pid].r/P[pid].r/2./P[pid].K));
}

double likelihood(double pars[])
{
  int i,j;
  double lnL=0.,p;
  int bincount[NBIN_R]={0},bincount_all[NBIN_R]={0};
//   decode_NFWprof(Z0,pow(10,pars[0])*M0,pow(10.,pars[1])*C0,VIR_C200,&Halo);  
  decode_NFWprof2(Z0,pow(10.,pars[0])*Rhos0,pow(10.,pars[1])*Rs0,VIR_C200,&Halo);

  #pragma omp parallel firstprivate(bincount) 
  {
    #pragma omp for
    for(i=0;i<nP;i++)
      solve_radial_orbit(i);
    
 #if ESTIMATOR==RADIAL_PHASE_ESTIMATOR
    #pragma omp for private(j)
    for(i=0;i<nP;i++)
    {  
      j=floor(NBIN_R*P[i].theta);
      if(j<0) j=0;
      if(j>=NBIN_R) j=NBIN_R-1;
      bincount[j]++;
    }
    #pragma omp critical  //note this is not barriered by default!!
    {
    for(j=0;j<NBIN_R;j++)
      bincount_all[j]+=bincount[j];
    }
    #pragma omp barrier  //this is necessary since OMP_critical has no implicit barrier!
// #pragma omp single private(j)
// {for(j=0;j<NBIN_R;j++) printf("%d,", bincount_all[j]);
// printf("\n");}
    #pragma omp for reduction(+:lnL)
    for(j=0;j<NBIN_R;j++)
    {
      lnL-=log_factorial(bincount_all[j]);
    }
#else
    #pragma omp for private(i,p,j,bincount) reduction(+:lnL)
    for(i=0;i<nP;i++)
    {     
  #if ESTIMATOR==MIXED_RADIAL_ESTIMATOR
      for(p=0,j=0;j<nP;j++)
      {
	if(P[i].r<P[j].rlim[0]||P[i].r>P[j].rlim[1]) continue;
	p+=vr_inv_part(P[i].r,j)/P[j].T;//contribution from j to i;
      }
      lnL+=log(p);
  #elif ESTIMATOR==ITERATIVE_ESTIMATOR
      lnL+=log(vr_inv_part(P[i].r,i)/P[i].T); 
  #elif ESTIMATOR==ENTROPY_ESTIMATOR
      lnL-=log(P[i].T);//maximum entropy estimator    
  #endif      
    }
#endif    
  }
//   printf("%g,%g: %g\n", pars[0],pars[1],lnL);
  return lnL;
}

int main(int argc, char **argv)
{
  double pars[NUM_PAR_MAX]={0.,0.};
  if(argc==3)
  {
    pars[0]=atof(argv[1]);
    pars[1]=atof(argv[2]);
  }
  
  init();
  freeze_energy(pars);
//  double x;
//  for(x=-1;x<1;x+=0.1)
//  {
//    pars[0]=x;
//    printf("%g,%g\n",pow(10.,x), likelihood(pars));
//  }
  double lnL=likelihood(pars);
  printf("Rhos=%g,Rs=%g,lnL=%g\n",pow(10.,pars[0]),pow(10.,pars[1]),lnL);

//   pars[0]+=0.1;
  freeze_and_like(pars);
  free_integration_space();
  
  return 0;
}

void init()
{
  char datafile[1024]=ROOTDIR"/data/mockhalo.hdf5";
  
  load_data(datafile);
  alloc_integration_space();
}
void freeze_energy(double pars[])
{//fix the energy parameter according to initial potential
  int i;
  decode_NFWprof2(Z0,pow(10.,pars[0])*Rhos0,pow(10.,pars[1])*Rs0,VIR_C200,&Halo);
  #pragma omp parallel for
  for(i=0;i<nP;i++)
      P[i].E=P[i].K+halo_pot(P[i].r);
}
double freeze_and_like(double pars[])
{
  freeze_energy(pars);
  return likelihood(pars);
}