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
	
#define EPS 1e-16
#define MODEL_MAX_INTVAL 1000
int NumericalIntegralRoutines=GSL_INTEGRAL;
double MODEL_TOL_BIN=1e-6, MODEL_TOL_BIN_ABS=1e-6, MODEL_TOL_REL=1e-3; //should be sufficient, good enough to constrain mass to 1% accuracy with 100000 particles
//accuracy in theta and phase-TS are approximately MODEL_TOL_REL, independent of nP.
//but it's still not accurate enough for minuit to work with the hessian; better use fmin(). TODO: implement MINUIT2, use SET EPSMachine to handle floating point precision
int RLimValConv=0; //convergence in value, not bin
int RLimInnerVal=0; //Adopt inner estimate for bin
// int NumericalIntegralRoutines=NR_INTEGRAL;
// double MODEL_TOL_BIN=1e-5, MODEL_TOL_BIN_ABS=1e-8, MODEL_TOL_REL=1e-6;

double HaloM0,HaloC0,HaloRhos0,HaloRs0,HaloZ0=0.;
struct NFWParZ Halo;

void choose_integral_routines(int type, int rlimvalconv, int rliminnerval)
{
  switch(type)
  {
	case NR_INTEGRAL:
	  MODEL_TOL_BIN=1e-6;
	  MODEL_TOL_BIN_ABS=1e-6;
	  MODEL_TOL_REL=1e-6;
	  break;
	case GSL_INTEGRAL:
	  MODEL_TOL_BIN=1e-6;//this is cheap
	  MODEL_TOL_BIN_ABS=1e-6;
	  MODEL_TOL_REL=1e-3;
	  break;
	default:
	  fprintf(stderr, "Error: unknown NumericalIntegralRoutines type=%d. must be 1 or 2\n", type);
	  exit(type);
  }
	NumericalIntegralRoutines=type;
	RLimValConv=rlimvalconv;
	RLimInnerVal=rliminnerval;
}
void define_halo(double pars[])
{ 
#if FIT_PAR_TYPE==PAR_TYPE_M_C
  decode_NFWprof(HaloZ0,pars[0]*HaloM0,pars[1]*HaloC0,VIR_C200,&Halo);  
#elif  FIT_PAR_TYPE==PAR_TYPE_RHOS_RS
  decode_NFWprof2(HaloZ0,pars[0]*Rhos0,pars[1]*HaloRs0,VIR_C200,&Halo);
#endif
}

double halo_pot(double r)
{
  double x=r/Halo.Rs;
  if(x<EPS) return Halo.Pots; //to avoid numerical nan at x=0;
  return Halo.Pots*log(1+x)/x; //the halo boundary is irrelevant in this problem, since only the potential difference affects the orbit
}
static int is_forbidden(double r, double E, double L2)
{//according to current potential and E,L,r
  if(-E-L2/2./r/r-halo_pot(r)<0) 
    return 1;

  return 0;
}

void solve_radial_limits ( Particle_t *P, double rmin, double rmax)
{
  //estimate pericenter and apocenter (outer bounds, to serve as integration limits),
  //according to current potential and E,L,r
  if ( is_forbidden ( P->r, P->E, P->L2 ) ) { //forbidden due to large difference between current pot and initial pot
	P->rlim[0]=P->r;
	P->rlim[1]=P->r;
	return;
  }
  
  char allow_left=0,allow_right=0;
  P->rlim[0]=rmin;
  P->rlim[1]=rmax;
  if ( !is_forbidden (P->rlim[0], P->E, P->L2 ) ) allow_left=1;
  if ( !is_forbidden (P->rlim[1], P->E, P->L2 ) ) allow_right=1;
  if (allow_left && allow_right) return; //no need to zoom
  
  double lbin[2], rbin[2], dx, xmid;
  if (!allow_left) //shrink left
  {
	lbin[0]=P->rlim[0];
	lbin[1]=P->r;
	dx=lbin[1]-lbin[0];
	while(1)
	{
	  dx/=2.;
	  xmid= ( lbin[0]+lbin[1] ) /2.;
	  if ( is_forbidden (xmid, P->E, P->L2) )
		lbin[0]=xmid;
	  else
		lbin[1]=xmid;
	  if(RLimValConv){
	  if (dx/lbin[0]<MODEL_TOL_BIN||dx<MODEL_TOL_BIN_ABS) break;
	  }else{
	  if (dx/( P->r-lbin[0])<MODEL_TOL_BIN||dx<MODEL_TOL_BIN_ABS) break;
	  }
	}
	if(RLimInnerVal)//no good, abaddon this
	P->rlim[0]=lbin[1];
	else
	P->rlim[0]=lbin[0];
  }
  if(!allow_right)//shrink right
  {
	rbin[0]=P->r;
	rbin[1]=P->rlim[1];
	dx=rbin[1]-rbin[0];
	while(1)
	{
	  dx/=2.;
	  xmid= ( rbin[0]+rbin[1] ) /2.;
	  if ( is_forbidden ( xmid, P->E, P->L2 ) )
		rbin[1]=xmid;
	  else
		rbin[0]=xmid;
	  if(RLimValConv){
	  if ( dx/rbin[1] <MODEL_TOL_BIN||dx<MODEL_TOL_BIN_ABS ) break;
	  }else
	  if ( dx/(rbin[1]-P->r ) <MODEL_TOL_BIN||dx<MODEL_TOL_BIN_ABS ) break;
	}
	if(RLimInnerVal)
	P->rlim[1]=rbin[0];
	else
	P->rlim[1]=rbin[1];
  }
}

double vr_inv_part(double r, double E, double L2)
{//E: binding energy, -(K+psi)
  double vr2=2*(-E-L2/2./r/r-halo_pot(r));
  if(vr2<=0) return 0.;
  return 1./sqrt(vr2);
}

typedef struct
{
  double E;
  double L2;
} OrbitPar;

static double vr_inv_rfunc(double r, void *params)
{ /*--- 1/vel_r ---*/
  return vr_inv_part(r, ((OrbitPar *)params)->E, ((OrbitPar *)params)->L2);
}
// static gsl_integration_workspace * GSL_workspace;
// #pragma omp threadprivate(GSL_workspace)
static gsl_integration_cquad_workspace * GSL_workspaceC;
#pragma omp threadprivate(GSL_workspaceC)
void alloc_integration_space()
{//have to allocate two workspaces when evaluating double integral, to avoid entangling inner and outer workspaces.
	#pragma omp parallel
  { //GSL_workspace=gsl_integration_workspace_alloc(MODEL_MAX_INTVAL);
	GSL_workspaceC=gsl_integration_cquad_workspace_alloc(MODEL_MAX_INTVAL);}
}
void free_integration_space()
{
	#pragma omp parallel
  { //gsl_integration_workspace_free (GSL_workspace);
    gsl_integration_cquad_workspace_free (GSL_workspaceC);}
}
void solve_radial_orbit(Particle_t *P, double rmin, double rmax, int estimator)//bottleneck in gsl_integration_cquad
{//find peri(apo)-centers and integrate the period 
  solve_radial_limits(P, rmin, rmax); //according to current potential and E,L,r; E must be initialized with a initial potential
  if(is_forbidden(P->r,P->E,P->L2))
  {
    P->T=1.;  //arbitrary, just to avoid 1/v/T=NaN, since 1/v=0.
    return;
  }
  
  gsl_function F;
  OrbitPar Fpar;
  Fpar.E=P->E;
  Fpar.L2=P->L2;
  F.function = &vr_inv_rfunc;
  F.params = &Fpar;
  double error;
  size_t neval;
  if(NumericalIntegralRoutines==NR_INTEGRAL)
  {
	double mid=P->rlim[0]+(P->rlim[1]-P->rlim[0])/2;
	P->T=qromo_gsl(&F, P->rlim[0], mid , midsql_gsl, MODEL_TOL_REL)+qromo_gsl(&F, mid, P->rlim[1], midsqu_gsl, MODEL_TOL_REL);
  }
  else
  {
	gsl_integration_cquad (&F, P->rlim[0], P->rlim[1], 0, MODEL_TOL_REL, //3, 
			  GSL_workspaceC, &(P->T), &error, &neval);
  }
  if(P->T<=0||isnan(P->T)||P->T==INFINITY){ fprintf(stderr,"Warning: T=%g, reset to 1. [%g,%g]\n", P->T, P->rlim[0], P->rlim[1]);P->T=1.;}
  if(IS_PHASE_ESTIMATOR(estimator))
  {
  double t;
  if(NumericalIntegralRoutines==NR_INTEGRAL)
  { 
	if(P->r<P->rlim[1])  
	  t=qromo_gsl(&F, P->rlim[0], P->r , midsql_gsl, MODEL_TOL_REL);
	else
	  t=P->T;
  }
  else
  {
	if(P->r<P->rlim[1])  
	gsl_integration_cquad (&F, P->rlim[0],P->r, MODEL_TOL_REL*P->T, MODEL_TOL_REL, //3, 
			  GSL_workspaceC, &t, &error, &neval);
	else
	  t=P->T;
  }
  #if PHASE_PERIOD==HALF_ORBIT_PERIOD
  P->theta=t/P->T; //AD test is sensitive to tails. this definition suits it best.
#elif PHASE_PERIOD==FULL_ORBIT_PERIOD //circular definition
  t=t/P->T/2.;
  P->theta=P->vr>=0?t:(1-t);//radial phase
#endif
	if(P->theta<=0.){
		//fprintf(stderr,"Warning: theta=%g (%g/%g), reset to EPS.\n", P->theta, t, P->T); 
		P->theta=EPS;}  //to avoid numerical problems
	if(P->theta>=1.){
		//fprintf(stderr,"Warning: theta=%g (%g,%g), reset to 1-EPS.\n", P->theta, t, P->T); 
		P->theta=1-EPS;}
//   if(P->theta==INFINITY) 
//   {printf("p=%d, r=%g, rlim=[%g,%g], t=%g, T=%g\n", pid, P->r, P->rlim[0], P->rlim[1], t, P->T*2);
//     printf("1/v=%g, isforbiidden=%d\n", vr_inv_part(P->rlim[0],pid), is_forbidden(P->rlim[0],pid));};
  }  
  //   result=smpintD(&F,xlim[0],xlim[1],MODEL_TOL_REL); //too slow
//   if(P->T<=0) printf("Part %d (M=%g,c=%g): r=%g,K=%g, E=%g, L2=%g; T=%g (vt/v=%g)\n",pid, Halo.M, Halo.c, P->r, P->K, P->E, P->L2, P->T, sqrt(P->L2/P->r/P->r/2./P->K));
}
double like_circular_moment(Tracer_t *Sample)
{//http://en.wikipedia.org/wiki/Circular_uniform_distribution
  //this resultant is rotational invariant
  int i,k=1;//k-th order moments, zoom in to examine the rotational symmetry on 2pi/k scale
  double c=0.,s=0.;
   #pragma omp parallel for reduction(+:c,s)
    for(i=0;i<Sample->nP;i++)
    {
      c+=cos(Sample->P[i].theta*2*k*M_PI);
      s+=sin(Sample->P[i].theta*2*k*M_PI);
    }
  return -2*k*(c*c+s*s)/Sample->nP;//chisquare(2) distributed; make negative to be comparable to likelihood
}
double like_cos_mean(Tracer_t *Sample)
{//http://en.wikipedia.org/wiki/Circular_uniform_distribution
  //this resultant is rotational invariant
  int i,k=1;//k-th order moments, zoom in to examine the rotational symmetry on 2pi/k scale
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<Sample->nP;i++)
      c+=cos(Sample->P[i].theta*2*k*M_PI);
  return -2*k*(c*c)/Sample->nP;//chisquare(2) distributed; make negative to be comparable to likelihood
}
double like_mean_phase(Tracer_t *Sample)
{  //raw mean phase (not squared)
  int i,k=1;//k-th order moments, zoom in to examine the rotational symmetry on 2pi/k scale
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<Sample->nP;i++)
    {
      c+=Sample->P[i].theta;
    }
    c/=Sample->nP;c-=0.5; 
    return c*sqrt(12.*Sample->nP);//standard normal variable
}
double like_linear_moment(Tracer_t *Sample)
{  
  int i,k=1;//k-th order moments, zoom in to examine the rotational symmetry on 2pi/k scale
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<Sample->nP;i++)
    {
      c+=Sample->P[i].theta;
    }
    c/=Sample->nP;c-=0.5; 
  //squared mean phase
    return -c*c*12*Sample->nP; //a standard chisquare variable
}

double KSTest(int FlagKuiper, Tracer_t *Sample)
{/*KS or Kuiper Test
  *KS test: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm
  *FlagKuiper=0: KS test
  *           1: Kuiper's test
  * return the probability of extreme values (p-value, P(TS>TSobs))
  */
  
  int i;
  double TS=0.,DL=0.,DU=0.,DLall=0.,DUall=0.,*theta;
  theta=malloc(sizeof(double)*Sample->nP);
  #pragma omp parallel firstprivate(DL,DU)
  {
    #pragma omp for 
    for(i=0;i<Sample->nP;i++)
      theta[i]=Sample->P[i].theta;
    #pragma omp single
    qsort(theta, Sample->nP, sizeof(double), cmpDouble);
    #pragma omp for
    for(i=0;i<Sample->nP;i++)
    {
      double a,b;
      a=theta[i]-i/(double)Sample->nP;
      b=(i+1.)/(double)Sample->nP-theta[i];
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
//   printf("%g: %g,%g\n", TS, KSprob(Sample->nP,TS), KSprobNR(Sample->nP,TS));
  double en=sqrt(Sample->nP);
  if(FlagKuiper)
#ifdef RETURN_PROB
    return KuiperProb(Sample->nP, TS);
#else
    return -(en+0.155+0.24/en)*TS;
#endif
  else
    #ifdef RETURN_PROB
    // return  lnL*sqrt(Sample->nP);
    return KSprob(Sample->nP, TS); //convert to something comparable to a log-likelihood
#else
    return -(en+0.12+0.11/en)*TS;
#endif
}
double AndersonDarlingTest_Old(Tracer_t *Sample)
{//Beloborodov&Levin 2004, apj, 613:224-237; deprecated. used the new function.
  int i;
  double lnL=0.,*theta;
  theta=malloc(sizeof(double)*Sample->nP);
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<Sample->nP;i++)
      theta[i]=Sample->P[i].theta;
    #pragma omp single
      {
      qsort(theta, Sample->nP, sizeof(double), cmpDouble);
      lnL=-log((1-theta[0])*theta[Sample->nP-1])-theta[0]+theta[Sample->nP-1]-1; //head and tail partitions
      }
    #pragma omp for reduction(+:lnL)
    for(i=1;i<Sample->nP;i++)
      lnL+=log(theta[i]/theta[i-1])*i*i/Sample->nP/Sample->nP
	  -log((1-theta[i])/(1-theta[i-1]))*(1.-(double)i/Sample->nP)*(1.-(double)i/Sample->nP)
	  +theta[i-1]-theta[i];
  }
  free(theta);
  lnL*=-Sample->nP; //differ by -1, to make it comparable to a loglike
  return lnL;
}
double AndersonDarlingTest(Tracer_t *Sample)
{//Beloborodov&Levin 2004, apj, 613:224-237; simplified equation as in the note.
  int i;
  double lnL=0., *theta;
  theta=malloc(sizeof(double)*Sample->nP);
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<Sample->nP;i++)
	{
      theta[i]=Sample->P[i].theta;
	}
    #pragma omp single
      qsort(theta, Sample->nP, sizeof(double), cmpDouble);
    #pragma omp for reduction(+:lnL)
    for(i=0;i<Sample->nP;i++)
      lnL+=-(1.+2.*i)*log(theta[i])+(1+2.*(i-Sample->nP))*log(1-theta[i]);
  }
  free(theta);
  lnL=lnL/Sample->nP-Sample->nP; 
  return -lnL; //differ by -1, to make it comparable to a loglike
}
double AndersonDarlingLike(Tracer_t *Sample, int estimator)
{/* adopting a PDF fit of ln(AD),
    return the log(Prob), the real likelihood (posterior probability)
    */
    double lnAD=log(-AndersonDarlingTest(Sample));
	static double gpar[2][3]={{0.569,   -0.570,    0.511}, { 0.431,    0.227,    0.569}}; //w, mu, sigma
	switch(estimator)
	{
	  case RADIAL_PHASE_AD_GEV:
		return log(GeneralizedExtremeValuePDF(lnAD, LnAD_GEV_MU, LnAD_GEV_SIGMA, LnAD_GEV_K)); //loglikelihood
	  case RADIAL_PHASE_AD_BINORMAL:
		return log(gpar[0][0]*NormPDF(lnAD, gpar[0][1], gpar[0][2])
				  +gpar[1][0]*NormPDF(lnAD, gpar[1][1], gpar[1][2]));
	  case RADIAL_PHASE_AD_NORMAL:
		lnAD=(lnAD-LnADMean)/LnADSig;
		return -lnAD*lnAD/2.; //loglike, subject to a const difference
	  default:
		DEBUGPRINT("Error: %d is not a ADLike estimator\n", estimator);
		exit(1);
	}
}
double like_phase_binned(Tracer_t *Sample)
{
  int i,j;
  double lnL=0.;
  int *bincount_all=calloc(Sample->nbin_r,sizeof(int));
 
#pragma omp parallel
  {
	int *bincount=calloc(Sample->nbin_r,sizeof(int));
    #pragma omp for private(j)
    for(i=0;i<Sample->nP;i++)
    {  
      j=floor(Sample->nbin_r*Sample->P[i].theta);
      if(j<0) j=0;
      if(j>=Sample->nbin_r) j=Sample->nbin_r-1;
      bincount[j]++;
    }
    #pragma omp critical  //note this is not barriered by default!!
    {
    for(j=0;j<Sample->nbin_r;j++)
      bincount_all[j]+=bincount[j];
    }
    #pragma omp barrier  //this is necessary since OMP_critical has no implicit barrier!
// #pragma omp single private(j)
// {for(j=0;j<Sample->nbin_r;j++) printf("%d,", bincount_all[j]); printf("\n");}
    #pragma omp for reduction(+:lnL)
    for(j=0;j<Sample->nbin_r;j++)
    {
      lnL-=log_factorial(bincount_all[j]);
    }
  }
  return lnL;
}
double like_phase_partition(Tracer_t *Sample)
{
  int i;
  double lnL=0.,*theta;
  theta=malloc(sizeof(double)*Sample->nP);
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<Sample->nP;i++)
      theta[i]=Sample->P[i].theta;
    #pragma omp single
      {
      qsort(theta, Sample->nP, sizeof(double), cmpDouble);
      lnL=log(theta[0]+1.-theta[Sample->nP-1]); //last partition
      }
    #pragma omp for reduction(+:lnL)
    for(i=1;i<Sample->nP;i++)
    {
      lnL+=log(theta[i]-theta[i-1]);
//       if(theta[i]==theta[i-1])
// 	printf("%d: %g,%g\n", i, theta[i], theta[i-1]);
    }
  }
  free(theta);
  return lnL;
}
double like_phase_process(Tracer_t *Sample)
{//junk
  int i;
  double lnL=0.;
#pragma omp parallel for reduction(+:lnL)
  for(i=0;i<Sample->nP;i++)
    lnL-=log(1-exp(Sample->nP*(Sample->P[i].theta-1)));
  
  return lnL;
}

void predict_radial_count (double RadialCountPred[], int nbin, Tracer_t *Sample)
{
    int i,j;
    double dr= ( Sample->rmax-Sample->rmin ) /nbin;

    for ( i=0; i<nbin; i++ ) {
        double rbin[2];
        rbin[0]=Sample->rmin+dr*i;
        rbin[1]=rbin[0]+dr;
        gsl_function F;
        F.function = &vr_inv_rfunc;
        double error,t;
        size_t neval;
        double p;
        p=0.;
        #pragma omp parallel for private(error,t,neval) firstprivate(F) reduction(+:p)
        for ( j=0; j<Sample->nP; j++ ) {
            if ( Sample->P[j].rlim[0]>=rbin[1]||Sample->P[j].rlim[1]<=rbin[0] ) continue;
			OrbitPar Fpar;
			Fpar.E=Sample->P[j].E;
			Fpar.L2=Sample->P[j].L2;
			F.params = &Fpar;
            gsl_integration_cquad ( &F, MAX ( rbin[0],Sample->P[j].rlim[0] ), MIN ( rbin[1],Sample->P[j].rlim[1] ), 0, MODEL_TOL_REL, //3,
                                    GSL_workspaceC, &t, &error, &neval );
            p+=t/Sample->P[j].T;
        }
        RadialCountPred[i]=p;
    }
}
double like_radial_bin(Tracer_t *Sample)
{
  int i;
  double lnL=0., dr=(Sample->rmax-Sample->rmin)/Sample->nbin_r;
  
    #pragma omp parallel for reduction(+:lnL)
    for(i=0;i<Sample->nbin_r;i++)
    { 
      if(0==Sample->RadialCount[i]) continue; //definitely skip empty bins (they have no contribution, but may produce NaNs)
      double p;
	  int j;
      //midpoint approximation
//       double r=dr*(i+0.5)+R_MIN;
//       for(p=0.,j=0;j<Sample->nP;j++)
//       {
// 	if(r<Sample->P[j].rlim[0]||r>Sample->P[j].rlim[1]) continue;
// 	p+=vr_inv_part(r,j)/Sample->P[j].T;//contribution from j to i;
//       }
//       p*=dr;
      //more proper way is to integrate inside bin
      double rbin[2]; rbin[0]=Sample->rmin+dr*i;rbin[1]=rbin[0]+dr;
      gsl_function F; F.function = &vr_inv_rfunc; 
	  OrbitPar Fpar;
      double error,t; size_t neval;
      for(p=0.,j=0;j<Sample->nP;j++)
      {
		if(Sample->P[j].rlim[0]>=rbin[1]||Sample->P[j].rlim[1]<=rbin[0]) continue;
		Fpar.E=Sample->P[j].E; 
		Fpar.L2=Sample->P[j].L2; 
		F.params = &Fpar;
		gsl_integration_cquad (&F, MAX(rbin[0],Sample->P[j].rlim[0]), MIN(rbin[1],Sample->P[j].rlim[1]), 0, MODEL_TOL_REL, //3, 
			       GSL_workspaceC, &t, &error, &neval);
		p+=t/Sample->P[j].T;
      }
      lnL+=Sample->RadialCount[i]*log(p);
    }
    
    return lnL;
}
double like_mixed_radial(Tracer_t *Sample)
{
  int i,j;
  double lnL=0.,p;
    #pragma omp parallel for private(p,j) reduction(+:lnL)
    for(i=0;i<Sample->nP;i++)
    {     
      for(p=0,j=0;j<Sample->nP;j++)
      {
		if(Sample->P[i].r<Sample->P[j].rlim[0]||Sample->P[i].r>Sample->P[j].rlim[1]) continue;
		p+=vr_inv_part(Sample->P[i].r,Sample->P[i].E, Sample->P[i].L2)/Sample->P[j].T;//contribution from j to i;
      }
      lnL+=log(p);
    }
    return lnL;
}
double like_entropy(Tracer_t *Sample)
{
    int i;
    double lnL=0.;
    #pragma omp parallel for reduction(+:lnL)
    for(i=0;i<Sample->nP;i++)
      lnL-=log(Sample->P[i].T);//maximum entropy estimator    
    return lnL;
}
double like_iterative_radial(Tracer_t *Sample)
{
    int i;
    double lnL=0.;
   #pragma omp parallel for reduction(+:lnL)
    for(i=0;i<Sample->nP;i++)
      lnL+=log(vr_inv_part(Sample->P[i].r, Sample->P[i].E, Sample->P[i].L2)/Sample->P[i].T); 
    return lnL;
}

void like_init(double pars[], int estimator, Tracer_t *Sample)
{
  int i;
  if(pars[0]<=0||pars[1]<=0) return;
  define_halo(pars);
  
  #pragma omp parallel for
  for(i=0;i<Sample->nP;i++)
    solve_radial_orbit(Sample->P+i,Sample->rmin,Sample->rmax,estimator);
}
double like_eval(double pars[], int estimator,Tracer_t *Sample)
{
  double lnL;
  
  if(pars[0]<=0||pars[1]<=0) return -INFINITY;
  switch(estimator)
  {
    case RADIAL_PHASE_BINNED:
      lnL=like_phase_binned(Sample);
      break;
    case RADIAL_PHASE_PARTITION:
      lnL=like_phase_partition(Sample);
      break;
    case RADIAL_PHASE_PROCESS:
      lnL=like_phase_process(Sample);
      break;
    case MIXED_RADIAL_ESTIMATOR:
      like_mixed_radial(Sample);
      break;
    case RADIAL_BIN_ESTIMATOR:
      lnL=like_radial_bin(Sample);
      break;    
    case ITERATIVE_ESTIMATOR:
      lnL=like_iterative_radial(Sample);
      break;
    case ENTROPY_ESTIMATOR:
      lnL=like_entropy(Sample);
      break;
    case RADIAL_PHASE_ROULETTE:
      lnL=AndersonDarlingTest(Sample);
      break;
	case RADIAL_PHASE_AD_GEV:
	case RADIAL_PHASE_AD_BINORMAL:
	case RADIAL_PHASE_AD_NORMAL:
      lnL=AndersonDarlingLike(Sample, estimator);
      break;  
    case RADIAL_PHASE_CMOMENT:
      lnL=like_circular_moment(Sample);
      break;
    case RADIAL_PHASE_COSMEAN:
      lnL=like_cos_mean(Sample);
      break;    
    case RADIAL_PHASE_LMOMENT:
      lnL=like_linear_moment(Sample);
      break;
	case RADIAL_PHASE_LMEANRAW:
      lnL=like_mean_phase(Sample);
      break;  
    case RADIAL_PHASE_KS:
      lnL=KSTest(0,Sample);  
      break;
    case RADIAL_PHASE_KUIPER:
      lnL=KSTest(1,Sample);
      break;
    default:
      DEBUGPRINT("Error: unknown Estimator=%d\n", estimator);
      exit(estimator);
  }  
  
  //   printf("%g,%g: %g\n", pars[0],pars[1],lnL);
  return lnL;
}

double likelihood(double pars[], int estimator, Tracer_t *Sample)
{
// 	time_t t1,t2,t3;
// 	t1=time(NULL);
  like_init(pars, estimator, Sample);
//   t2=time(NULL);
  double lnL=like_eval(pars, estimator, Sample);
//   t3=time(NULL);
//   printf("Time: %ld,%ld, %ld,%ld,%ld\n", t2-t1, t3-t2, t1,t2,t3);
  return lnL;
}

double like_to_chi2(double lnL, int estimator)
{//convert likelihood() values to a chi-square measure
  switch(estimator)
  {
	case RADIAL_PHASE_LMEANRAW:
	  return lnL*lnL;
	case RADIAL_PHASE_LMOMENT:
	  return -lnL;
	case RADIAL_PHASE_ROULETTE:
	  lnL=(log(-lnL)-LnADMean)/LnADSig;
	  return lnL*lnL;
	case RADIAL_PHASE_AD_GEV:
	case RADIAL_PHASE_AD_BINORMAL:
	case RADIAL_PHASE_AD_NORMAL:  
	case MIXED_RADIAL_ESTIMATOR:
	  return -2.*lnL; //simply 2*(the negative loglike), to make it comparable to chisquare; Note these could differ from the chi-2 by a constant.
// 	case RADIAL_BIN_ESTIMATOR:
// 	  return -lnL*2;  //simply 2*(the negative loglike), to make it comparable to chisquare
	default:
	  fprintf(stderr, "Error: like_to_chi2() not implemented for estimator=%d yet\n", estimator);
	  exit(1);
  }
}

void freeze_energy(double pars[], Tracer_t *Sample)
{//fix the energy parameter according to initial potential
  int i;
  define_halo(pars);
  #pragma omp parallel for
  for(i=0;i<Sample->nP;i++)
      Sample->P[i].E=-(Sample->P[i].K+halo_pot(Sample->P[i].r));//differ from previous version
}

double freeze_and_like(double pars[], int estimator, Tracer_t *Sample)
{
  freeze_energy(pars, Sample);
  return likelihood(pars, estimator, Sample);
}

double jointLE_Flike(double pars[], int estimator, int nbinL, int nbinE, Tracer_t *Sample)
{//automatically update views if needed; then freeze and like
  int i;
  double chi2;
  if(Sample->nView!=nbinL||Sample->ViewType!='L') //already allocated
	create_tracer_views(Sample, nbinL, 'L');
  for(i=0,chi2=0;i<nbinL;i++)
	chi2+=jointE_Flike(pars, estimator, nbinE, Sample->Views+i);
  return chi2;
}
double jointE_Flike(double pars[], int estimator, int nbin, Tracer_t *Sample)
{//this does freeze_and_like
  int i;
  double lnL, chi2;
  freeze_energy(pars, Sample);
  create_tracer_views(Sample, nbin, 'E');
  for(i=0,chi2=0;i<nbin;i++)
  {
// 	if(estimator==RADIAL_BIN_ESTIMATOR)  count_tracer_radial(View+i, NumRadialBin);
	lnL=likelihood(pars, estimator, Sample->Views+i);
	chi2+=like_to_chi2(lnL, estimator);
// 	if(estimator==RADIAL_BIN_ESTIMATOR)  free_tracer_rcount(View+i);
  }
  free_tracer_views(Sample);
  return chi2;
}
void create_nested_views(double pars[], int nbin[], char ViewTypes[], Tracer_t *Sample)//iteratively: freeze, fit; then freeze, then fit.
{//ViewTypes should be a non-empty string!
  //used to create static views (jointE, jointLE dynamically create views by themselves instead)
  //do not mix this with joint_like functions, since they destroy the views!
  int i;
  if(ViewTypes[0]=='E')  freeze_energy(pars, Sample);
  create_tracer_views(Sample, nbin[0], ViewTypes[0]);
  if(ViewTypes[1]=='\0') return; //done
  
  for(i=0;i<nbin[0];i++)
	create_nested_views(pars, nbin+1, ViewTypes+1, Sample->Views+i);//descend
}
double nested_views_like(double pars[], int estimator, Tracer_t *Sample)
{//pure like, without freezing energy
  double lnL;
  if(!Sample->nView)
  {
	lnL=likelihood(pars, estimator, Sample);
	return like_to_chi2(lnL,estimator);
  }
  
  int i;
  for(i=0,lnL=0.;i<Sample->nView;i++)
	lnL+=nested_views_like(pars, estimator, Sample->Views+i);
  return lnL;
}
double nested_views_Flike(double pars[], int estimator, Tracer_t *Sample)
{//freeze and like
  freeze_energy(pars, Sample);
  return nested_views_like(pars, estimator, Sample);
}