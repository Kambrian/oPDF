/*ToDo:
 * Same data for good estimators, directly compare their performance (likelihood curve, CI) (3.84/2 for log-like, -2.48 for roulette, at 95% CL)
 * shift to new data
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mymath.h"
#include "cosmology.h"
#include "io.h"
#include "models.h"

#define Z0 0. 
double M0,C0,Rhos0,Rs0;
struct NFWParZ Halo;

double halo_pot(double r)
{
  double x=r/Halo.Rs;
  if(x<1e-16) return Halo.Pots; //to avoid numerical nan at x=0;
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
  if(vr2<=0) return 0.;
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
void solve_radial_orbit(int pid, int estimator)
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
  if(IS_PHASE_ESTIMATOR(estimator))
  {
  double t;
  gsl_integration_cquad (&F, P[pid].rlim[0],P[pid].r, 0, MODEL_TOL_REL, //3, 
			 GSL_workspaceC, &t, &error, &neval);
#if PHASE_PERIOD==HALF_ORBIT_PERIOD
  P[pid].theta=t/P[pid].T; //AD test is sensitive to tails. this definition suits it best.
#elif PHASE_PERIOD==FULL_ORBIT_PERIOD //circular definition
  t=t/P[pid].T/2.;
  P[pid].theta=P[pid].vr>=0?t:(1-t);//radial phase
#endif
//   if(P[pid].theta==INFINITY) 
//   {printf("p=%d, r=%g, rlim=[%g,%g], t=%g, T=%g\n", pid, P[pid].r, P[pid].rlim[0], P[pid].rlim[1], t, P[pid].T*2);
//     printf("1/v=%g, isforbiidden=%d\n", vr_inv_part(P[pid].rlim[0],pid), is_forbidden(P[pid].rlim[0],pid));};
  }  
  //   result=smpintD(&F,xlim[0],xlim[1],MODEL_TOL_REL); //too slow
//   if(P[pid].T<=0) printf("Part %d (M=%g,c=%g): r=%g,K=%g, E=%g, L2=%g; T=%g (vt/v=%g)\n",pid, Halo.M, Halo.c, P[pid].r, P[pid].K, P[pid].E, P[pid].L2, P[pid].T, sqrt(P[pid].L2/P[pid].r/P[pid].r/2./P[pid].K));
}
double like_circular_moment()
{//http://en.wikipedia.org/wiki/Circular_uniform_distribution
  //this resultant is rotational invariant
  int i,k=1;//k-th order moments, zoom in to examine the rotational symmetry on 2pi/k scale
  double c=0.,s=0.;
   #pragma omp parallel for reduction(+:c,s)
    for(i=0;i<nP;i++)
    {
      c+=cos(P[i].theta*2*k*M_PI);
      s+=sin(P[i].theta*2*k*M_PI);
    }
  return -2*k*(c*c+s*s)/nP;//chisquare(2) distributed; make negative to be comparable to likelihood
}
double like_cos_mean()
{//http://en.wikipedia.org/wiki/Circular_uniform_distribution
  //this resultant is rotational invariant
  int i,k=1;//k-th order moments, zoom in to examine the rotational symmetry on 2pi/k scale
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<nP;i++)
      c+=cos(P[i].theta*2*k*M_PI);
  return -2*k*(c*c)/nP;//chisquare(2) distributed; make negative to be comparable to likelihood
}
double like_linear_moment()
{  
  int i,k=1;//k-th order moments, zoom in to examine the rotational symmetry on 2pi/k scale
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<nP;i++)
    {
      c+=P[i].theta;
    }
    c/=nP;c-=0.5; 
#ifdef RETURN_RAWMEAN
    return c*sqrt(12.*nP);//standard normal variable
#else    
    return -c*c*12*nP; //a standard chisquare variable
#endif
}

double KSTest(int FlagKuiper)
{/*KS or Kuiper Test
  *KS test: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm
  *FlagKuiper=0: KS test
  *           1: Kuiper's test
  * return the probability of extreme values (p-value, P(TS>TSobs))
  */
  
  int i;
  double TS=0.,DL=0.,DU=0.,DLall=0.,DUall=0.,*theta;
  theta=malloc(sizeof(double)*nP);
  #pragma omp parallel firstprivate(DL,DU)
  {
    #pragma omp for 
    for(i=0;i<nP;i++)
      theta[i]=P[i].theta;
    #pragma omp single
    qsort(theta, nP, sizeof(double), cmpDouble);
    #pragma omp for
    for(i=0;i<nP;i++)
    {
      double a,b;
      a=theta[i]-i/(double)nP;
      b=(i+1.)/(double)nP-theta[i];
      if(DU<a) DU=a;
      if(DL<b) DL=b;
    }
    #pragma omp critical
    {
      if(DUall<DU) DUall=DU;
      if(DLall<DL) DLall=DL;
    }
    #pragma omp barrier
    #pragma omp single
    {
      if(FlagKuiper)
	TS=DLall+DUall;
      else
	TS=MAX(DUall,DLall);
    }
  }
  free(theta);
//   printf("%g: %g,%g\n", TS, KSprob(nP,TS), KSprobNR(nP,TS));
  double en=sqrt(nP);
  if(FlagKuiper)
#ifdef RETURN_PROB
    return KuiperProb(nP, TS);
#else
    return -(en+0.155+0.24/en)*TS;
#endif
  else
    #ifdef RETURN_PROB
    // return  lnL*sqrt(nP);
    return KSprob(nP, TS); //convert to something comparable to a log-likelihood
#else
    return -(en+0.12+0.11/en)*TS;
#endif
}
double AndersonDarlingTest_Old()
{//Beloborodov&Levin 2004, apj, 613:224-237; deprecated. used the new function.
  int i;
  double lnL=0.,*theta;
  theta=malloc(sizeof(double)*nP);
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<nP;i++)
      theta[i]=P[i].theta;
    #pragma omp single
      {
      qsort(theta, nP, sizeof(double), cmpDouble);
      lnL=-log((1-theta[0])*theta[nP-1])-theta[0]+theta[nP-1]-1; //head and tail partitions
      }
    #pragma omp for reduction(+:lnL)
    for(i=1;i<nP;i++)
      lnL+=log(theta[i]/theta[i-1])*i*i/nP/nP
	  -log((1-theta[i])/(1-theta[i-1]))*(1.-(double)i/nP)*(1.-(double)i/nP)
	  +theta[i-1]-theta[i];
  }
  free(theta);
  lnL*=-nP; //differ by -1, to make it comparable to a loglike
  return lnL;
}
double AndersonDarlingTest()
{//Beloborodov&Levin 2004, apj, 613:224-237; simplified equation as in the note.
  int i;
  double lnL=0., *theta;
  theta=malloc(sizeof(double)*nP);
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<nP;i++)
      theta[i]=P[i].theta;
    #pragma omp single
      qsort(theta, nP, sizeof(double), cmpDouble);
    #pragma omp for reduction(+:lnL)
    for(i=0;i<nP;i++)
      lnL+=-(1.+2.*i)*log(theta[i])+(1+2.*(i-nP))*log(1-theta[i]);
  }
  free(theta);
  lnL=lnL/nP-nP; 
  return -lnL; //differ by -1, to make it comparable to a loglike
}
double like_phase_binned()
{
  int i,j;
  double lnL=0.;
  int bincount[NBIN_R]={0},bincount_all[NBIN_R]={0};
  
#pragma omp parallel firstprivate(bincount)
  {
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
// {for(j=0;j<NBIN_R;j++) printf("%d,", bincount_all[j]); printf("\n");}
    #pragma omp for reduction(+:lnL)
    for(j=0;j<NBIN_R;j++)
    {
      lnL-=log_factorial(bincount_all[j]);
    }
  }
  return lnL;
}
double like_phase_partition()
{
  int i;
  double lnL=0.,*theta;
  theta=malloc(sizeof(double)*nP);
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<nP;i++)
      theta[i]=P[i].theta;
    #pragma omp single
      {
      qsort(theta, nP, sizeof(double), cmpDouble);
      lnL=log(theta[0]+1.-theta[nP-1]); //last partition
      }
    #pragma omp for reduction(+:lnL)
    for(i=1;i<nP;i++)
    {
      lnL+=log(theta[i]-theta[i-1]);
//       if(theta[i]==theta[i-1])
// 	printf("%d: %g,%g\n", i, theta[i], theta[i-1]);
    }
  }
  free(theta);
  return lnL;
}
double like_phase_process()
{//junk
  int i;
  double lnL=0.;
#pragma omp parallel for reduction(+:lnL)
  for(i=0;i<nP;i++)
    lnL-=log(1-exp(nP*(P[i].theta-1)));
  
  return lnL;
}

static int RadialCountAll[NBIN_R]={0};
void fill_radial_bin()
{//count inside linear radial bins
  int i,j;
  double dr=(R_MAX-R_MIN)/NBIN_R;
  int bincount[NBIN_R]={0};
  
  for(i=0;i<NBIN_R;i++)  RadialCountAll[i]=0; 
#pragma omp parallel firstprivate(bincount)
  {
    #pragma omp for private(j)
    for(i=0;i<nP;i++)
    {  
      j=floor((P[i].r-R_MIN)/dr);
      if(j<0) j=0;
      if(j>=NBIN_R) j=NBIN_R-1;
      bincount[j]++;
    }
    #pragma omp critical  //note this is not barriered by default!!
    {
    for(j=0;j<NBIN_R;j++)
      RadialCountAll[j]+=bincount[j];
    }
  }
}
double like_radial_bin()
{
  int i,j;
  double lnL=0., dr=(R_MAX-R_MIN)/NBIN_R;
  
    #pragma omp parallel for private(j) reduction(+:lnL)
    for(i=0;i<NBIN_R;i++)
    { 
      if(0==RadialCountAll[i]) continue; //definitely skip empty bins (they have no contribution, but may produce NaNs)
      double p;
      //midpoint approximation
//       double r=dr*(i+0.5)+R_MIN;
//       for(p=0.,j=0;j<nP;j++)
//       {
// 	if(r<P[j].rlim[0]||r>P[j].rlim[1]) continue;
// 	p+=vr_inv_part(r,j)/P[j].T;//contribution from j to i;
//       }
//       p*=dr;
      //more proper way is to integrate inside bin
      double rbin[2]; rbin[0]=R_MIN+dr*i;rbin[1]=rbin[0]+dr;
      gsl_function F; F.function = &vr_inv_rfunc; F.params = &j;
      double error,t; size_t neval;
      for(p=0.,j=0;j<nP;j++)
      {
	if(P[j].rlim[0]>=rbin[1]||P[j].rlim[1]<=rbin[0]) continue;
	gsl_integration_cquad (&F, MAX(rbin[0],P[j].rlim[0]), MIN(rbin[1],P[j].rlim[1]), 0, MODEL_TOL_REL, //3, 
			       GSL_workspaceC, &t, &error, &neval);
	p+=t/P[j].T;
      }
      lnL+=RadialCountAll[i]*log(p);
    }
    
    return lnL;
}
double like_mixed_radial()
{
  int i,j;
  double lnL=0.,p;
    #pragma omp parallel for private(p,j) reduction(+:lnL)
    for(i=0;i<nP;i++)
    {     
      for(p=0,j=0;j<nP;j++)
      {
	if(P[i].r<P[j].rlim[0]||P[i].r>P[j].rlim[1]) continue;
	p+=vr_inv_part(P[i].r,j)/P[j].T;//contribution from j to i;
      }
      lnL+=log(p);
    }
    return lnL;
}
double like_entropy()
{
    int i;
    double lnL=0.;
    #pragma omp parallel for reduction(+:lnL)
    for(i=0;i<nP;i++)
      lnL-=log(P[i].T);//maximum entropy estimator    
    return lnL;
}
double like_iterative_radial()
{
    int i;
    double lnL=0.;
   #pragma omp parallel for reduction(+:lnL)
    for(i=0;i<nP;i++)
      lnL+=log(vr_inv_part(P[i].r,i)/P[i].T); 
    return lnL;
}

void like_init(double pars[], int estimator)
{
  int i;
  if(pars[0]<0||pars[1]<0) return;
  define_halo(pars);
  
  #pragma omp parallel for
  for(i=0;i<nP;i++)
    solve_radial_orbit(i,estimator);
}
double like_eval(double pars[], int estimator)
{
  double lnL;
  
  if(pars[0]<0||pars[1]<0) return -INFINITY;
  switch(estimator)
  {
    case RADIAL_PHASE_BINNED:
      lnL=like_phase_binned();
      break;
    case RADIAL_PHASE_PARTITION:
      lnL=like_phase_partition();
      break;
    case RADIAL_PHASE_PROCESS:
      lnL=like_phase_process();
      break;
    case MIXED_RADIAL_ESTIMATOR:
      like_mixed_radial();
      break;
    case RADIAL_BIN_ESTIMATOR:
      lnL=like_radial_bin();
      break;    
    case ITERATIVE_ESTIMATOR:
      lnL=like_iterative_radial();
      break;
    case ENTROPY_ESTIMATOR:
      lnL=like_entropy();
      break;
    case RADIAL_PHASE_ROULETTE:
      lnL=AndersonDarlingTest();
      break;
    case RADIAL_PHASE_CMOMENT:
      lnL=like_circular_moment();
      break;
    case RADIAL_PHASE_COSMEAN:
      lnL=like_cos_mean();
      break;    
    case RADIAL_PHASE_LMOMENT:
      lnL=like_linear_moment();
      break;
    case RADIAL_PHASE_KS:
      lnL=KSTest(0);  
      break;
    case RADIAL_PHASE_KUIPER:
      lnL=KSTest(1);
      break;
    default:
      fprintf(stderr, "Error: unknown Estimator=%d\n", estimator);
      exit(estimator);
  }  
  
  //   printf("%g,%g: %g\n", pars[0],pars[1],lnL);
  return lnL;
}

double likelihood(double pars[], int estimator)
{
  like_init(pars, estimator);
  return like_eval(pars, estimator);
}

void init()
{
  char datafile[1024]=ROOTDIR"/data/mockhalo_wenting.hdf5";
  SUBSAMPLE_SIZE=1000;
  R_MIN=1;
  R_MAX=1000;
  M0=183.5017;
  C0=16.1560;
  
  if(NULL!=getenv("DynDataFile"))
  {
    printf("Importing parameters from environment..\n");
    sprintf(datafile,"%s/data/%s", ROOTDIR, getenv("DynDataFile"));
    SUBSAMPLE_SIZE=strtol(getenv("DynSIZE"),NULL, 10);
    R_MIN=strtod(getenv("DynRMIN"),NULL);
    R_MAX=strtod(getenv("DynRMAX"),NULL);
    M0=strtod(getenv("DynM0"),NULL);
    C0=strtod(getenv("DynC0"),NULL);
  }
  else
    printf("Warning: Using default parameters with datafile %s .\n", datafile);
 
  printf("%s; %g,%g;%g,%g\n", datafile, R_MIN,R_MAX,M0,C0);
  decode_NFWprof(Z0,M0,C0,VIR_C200,&Halo);
  Rhos0=Halo.Rhos;
  Rs0=Halo.Rs;
  
  load_data(datafile);
  shuffle_data(100);
  alloc_integration_space();
}
void select_particles(int subsample_id)
{
  sample_data(subsample_id);
  fill_radial_bin();
}
void freeze_energy(double pars[])
{//fix the energy parameter according to initial potential
  int i;
  define_halo(pars);
  #pragma omp parallel for
  for(i=0;i<nP;i++)
      P[i].E=P[i].K+halo_pot(P[i].r);
}
void define_halo(double pars[])
{ 
#if FIT_PAR_TYPE==PAR_TYPE_M_C
  decode_NFWprof(Z0,pars[0]*M0,pars[1]*C0,VIR_C200,&Halo);  
#elif  FIT_PAR_TYPE==PAR_TYPE_RHOS_RS
  decode_NFWprof2(Z0,pars[0]*Rhos0,pars[1]*Rs0,VIR_C200,&Halo);
#endif
}
double freeze_and_like(double pars[], int estimator)
{
  freeze_energy(pars);
  return likelihood(pars, estimator);
}