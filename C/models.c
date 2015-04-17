#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>
#include <time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mymath.h"
#include "globals.h"
#include "cosmology.h"
#include "tracer.h"
#include "models.h"
	
#define MODEL_MAX_INTVAL 1000

static int is_forbidden(double r, double E, double L2, Halo_t *halo)
{//according to current potential and E,L,r
  if(-E-L2/2./r/r-halo_pot(r, halo)<0) 
    return 1;

  return 0;
}

void solve_radial_limits ( Particle_t *P, double rmin, double rmax, Halo_t *halo)
{
  //estimate pericenter and apocenter (outer bounds, to serve as integration limits),
  //according to current potential and E,L,r
  if ( is_forbidden ( P->r, P->E, P->L2, halo ) ) { //forbidden due to large difference between current pot and initial pot
	P->rlim[0]=P->r;
	P->rlim[1]=P->r;
	return;
  }
  
  char allow_left=0,allow_right=0;
  P->rlim[0]=rmin;
  P->rlim[1]=rmax;
  if ( !is_forbidden (P->rlim[0], P->E, P->L2, halo ) ) allow_left=1;
  if ( !is_forbidden (P->rlim[1], P->E, P->L2, halo ) ) allow_right=1;
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
	  if ( is_forbidden (xmid, P->E, P->L2, halo) )
		lbin[0]=xmid;
	  else
		lbin[1]=xmid;
	  if (dx/( P->r-lbin[0])<Globals.tol.bin||dx<Globals.tol.bin_abs) break;
	}
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
	  if ( is_forbidden ( xmid, P->E, P->L2, halo ) )
		rbin[1]=xmid;
	  else
		rbin[0]=xmid;
	  if ( dx/(rbin[1]-P->r ) <Globals.tol.bin||dx<Globals.tol.bin_abs ) break;
	}
	P->rlim[1]=rbin[1];
  }
}

double vr_inv_part(double r, double E, double L2, Halo_t *halo)
{//E: binding energy, -(K+psi)
  double vr2=2*(-E-L2/2./r/r-halo_pot(r, halo));//what about hubble flow? this is physical velocity. after excluding hubble flow, the comoving velocity would be asymmetric for in-out motions...
  if(vr2<=0) return 0.;
  return 1./sqrt(vr2);
}

typedef struct
{
  double E;
  double L2;
  Halo_t *halo;
} OrbitPar;

static double vr_inv_rfunc(double r, void *params)
{ /*--- 1/vel_r ---*/
  return vr_inv_part(r, ((OrbitPar *)params)->E, ((OrbitPar *)params)->L2, ((OrbitPar *)params)->halo);
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
void solve_radial_orbit(Particle_t *P, double rmin, double rmax, Halo_t *halo, int FlagSetPhase)//bottleneck in gsl_integration_cquad
{//find peri(apo)-centers and integrate the period 
  
  if(halo->IsForbidden)
	{
	  printf("Halo not allowed with param [%g,%g].\n Not solving orbits.\n", halo->pars[0], halo->pars[1]);
	  return;
	}
  solve_radial_limits(P, rmin, rmax, halo); //according to current potential and E,L,r; E must be initialized with a initial potential
  if(is_forbidden(P->r,P->E,P->L2, halo))
  {
    P->T=1.;  //arbitrary, just to avoid 1/v/T=NaN, since 1/v=0.
    return;
  }
  
  gsl_function F;
  OrbitPar Fpar;
  Fpar.E=P->E;
  Fpar.L2=P->L2;
  Fpar.halo=halo;
  F.function = &vr_inv_rfunc;
  F.params = &Fpar;
  
  double error;
  //   gsl_integration_qags (&F, xlim[0],xlim[1], 0, Globals.tol.rel, MODEL_MAX_INTVAL, //3, 
  // 		       GSL_workspace, &result, &error);
  size_t neval;
  gsl_integration_cquad (&F, P->rlim[0], P->rlim[1], 0, Globals.tol.rel, //3, 
			 GSL_workspaceC, &(P->T), &error, &neval);
  if(P->T<=0||isnan(P->T)||P->T==INFINITY){ 
	fprintf(stderr,"Warning: T=%g, reset to 1. [%g,%g]\n", P->T, P->rlim[0], P->rlim[1]);P->T=1.;}
  if(FlagSetPhase)
  {
  double t;
  gsl_integration_cquad (&F, P->rlim[0],P->r, 0, Globals.tol.rel, //3, 
			 GSL_workspaceC, &t, &error, &neval);
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
  //   result=smpintD(&F,xlim[0],xlim[1],Globals.tol.rel); //too slow
//   if(P->T<=0) printf("Part %d (M=%g,c=%g): r=%g,K=%g, E=%g, L2=%g; T=%g (vt/v=%g)\n",pid, Halo.M, Halo.c, P->r, P->K, P->E, P->L2, P->T, sqrt(P->L2/P->r/P->r/2./P->K));
}
double MeanPhaseTest(Tracer_t *Sample, int FlagRaw)
{ 
  if(Sample->halo.IsForbidden) return INFINITY;

  int i;
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<Sample->nP;i++)
    {
      c+=Sample->P[i].theta;
    }
    c/=Sample->nP;c-=0.5; 
	if(FlagRaw)  
	  return c*sqrt(12.*Sample->nP);//raw mean phase, standard normal variable
	else//squared mean phase
	  return c*c*12*Sample->nP; //a standard chisquare variable	  
}
double AndersonDarlingTest(Tracer_t *Sample)
{//Beloborodov&Levin 2004, apj, 613:224-237; simplified equation as in the note.
  if(Sample->halo.IsForbidden) return INFINITY;

  int i;
  double AD=0., *theta;
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
    #pragma omp for reduction(+:AD)
    for(i=0;i<Sample->nP;i++)
      AD+=-(1.+2.*i)*log(theta[i])+(1+2.*(i-Sample->nP))*log(1-theta[i]);
  }
  free(theta);
  AD=AD/Sample->nP-Sample->nP; 
  return AD; //AD test stat
}
void predict_radial_count (double RadialCountPred[], int nbin, int FlagRLogBin,  Tracer_t *Sample)
{
    if(Sample->halo.IsForbidden)
	{
	  printf("Halo not allowed with param [%g,%g]\n Not predicting counts.\n", Sample->halo.pars[0], Sample->halo.pars[1]);
	  return;
	}
	
    int i,j;
	double dr, logRmin, factor;
// 	printf("Log=%d\n",FlagRLogBin);
	if(FlagRLogBin)
	{
	  logRmin=log(Sample->rmin);
	  dr=(log(Sample->rmax)-logRmin)/nbin;
	  factor=exp(dr);
	}
	else
	  dr=(Sample->rmax-Sample->rmin)/nbin;

    for ( i=0; i<nbin; i++ ) {
        double rbin[2];
		if(FlagRLogBin)
	  {
		rbin[0]=exp(logRmin+dr*i);
		rbin[1]=rbin[0]*factor;
	  }
	  else
	  {
        rbin[0]=Sample->rmin+dr*i;
        rbin[1]=rbin[0]+dr;
	  }
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
			Fpar.halo=&(Sample->halo);
			F.params = &Fpar;
            gsl_integration_cquad ( &F, MAX ( rbin[0],Sample->P[j].rlim[0] ), MIN ( rbin[1],Sample->P[j].rlim[1] ), 0, Globals.tol.rel, //3,
                                    GSL_workspaceC, &t, &error, &neval );
            p+=t/Sample->P[j].T;
        }
        RadialCountPred[i]=p;
    }
}
double like_radial_bin( Tracer_t *Sample)
{
  if(Sample->halo.IsForbidden) return -INFINITY;

  int i;
  double lnL=0., dr, logRmin, factor;
  if(Sample->nbin_r==0) 
  {
	DEBUGPRINT("Error: %d radial bins. call count_tracer_radial() before radial like.\n", 0);
	exit(1);
  }
  if(Sample->FlagRLogBin)
  {
	logRmin=log(Sample->rmin);
	dr=(log(Sample->rmax)-logRmin)/Sample->nbin_r;
	factor=exp(dr);
  }
  else
	dr=(Sample->rmax-Sample->rmin)/Sample->nbin_r;
  
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
      double rbin[2];
	  if(Sample->FlagRLogBin)
	  {
		rbin[0]=exp(logRmin+dr*i);
		rbin[1]=rbin[0]*factor;
	  }
	  else
	  {
		rbin[0]=Sample->rmin+dr*i;
		rbin[1]=rbin[0]+dr;
	  }
      gsl_function F; F.function = &vr_inv_rfunc; 
	  OrbitPar Fpar;
      double error,t; size_t neval;
      for(p=0.,j=0;j<Sample->nP;j++)
      {
		if(Sample->P[j].rlim[0]>=rbin[1]||Sample->P[j].rlim[1]<=rbin[0]) continue;
		Fpar.E=Sample->P[j].E; 
		Fpar.L2=Sample->P[j].L2; 
		Fpar.halo=&(Sample->halo);
		F.params = &Fpar;
		gsl_integration_cquad (&F, MAX(rbin[0],Sample->P[j].rlim[0]), MIN(rbin[1],Sample->P[j].rlim[1]), 0, Globals.tol.rel, //3, 
			       GSL_workspaceC, &t, &error, &neval);
		p+=t/Sample->P[j].T;
      }
      lnL+=Sample->RadialCount[i]*log(p);
    }
    return lnL;
}

double like_eval(Estimator_t estimator, Tracer_t *Sample)
{
  double lnL;
  switch(estimator)
  {
    case EID_RBinLike:
      lnL=like_radial_bin(Sample);
      break;    
    case EID_PhaseAD:
      lnL=AndersonDarlingTest(Sample);
      break;
    case EID_PhaseMean:
	  lnL=MeanPhaseTest(Sample, 0);
	  break;
	case EID_PhaseMeanRaw:
      lnL=MeanPhaseTest(Sample, 1);
      break;
    default:
      DEBUGPRINT("Error: unknown Estimator=%d\n", estimator);
      exit(estimator);
  }
//   printf("[%d] %g,%g: %g\n", Sample->halo.IsForbidden, Sample->halo.pars[0],Sample->halo.pars[1],lnL);
  Sample->lnL=lnL;
  return lnL;
}

double nested_views_like(Estimator_t estimator,  Tracer_t *Sample, int nbin_r, int FlagRLogBin)
{//pure like, without freeze_energy()/set_orbits() (except for set_orbits() for r-Views);  prepare them before calling this!
  //nbin_r and FlagRLogBin are only used when estimator is RBinLike  
  double lnL;
  if(!Sample->nView)
  {
	if(estimator==EID_RBinLike)  count_tracer_radial(Sample, nbin_r, FlagRLogBin);
	lnL=like_eval(estimator, Sample);
	if(estimator==EID_RBinLike)  free_tracer_rcounts(Sample);
	return lnL;
  }
  
  int i;
  for(i=0,lnL=0.;i<Sample->nView;i++)
  {
	if(Sample->ViewType=='r')
	  tracer_set_orbits(Sample->Views+i, estimator!=EID_RBinLike);//update radial bounds
	lnL+=nested_views_like(estimator, Sample->Views+i, nbin_r, FlagRLogBin);
  }
  Sample->lnL=lnL;
  return lnL;
}

struct like_args
{
  Tracer_t *Sample;
  Estimator_t estimator;
};
static double like_rfunc(const gsl_vector *x, void *params)
{
  struct like_args *arg;
  arg=(struct like_args *)params;
  double pars[NUM_PAR_MAX];
  int i;
  for(i=0;i<x->size;i++)
	pars[i]=gsl_vector_get(x, i);
  halo_set_param(pars, &(arg->Sample->halo));
  tracer_set_energy(arg->Sample);
  if(arg->estimator==EID_RBinLike)
	tracer_set_orbits(arg->Sample, 0);
  else
	tracer_set_orbits(arg->Sample, 1);
  double lnL;
  lnL=like_eval(arg->estimator, arg->Sample);
  if(arg->estimator==EID_RBinLike)
	return -lnL;//to be minimized
  return lnL;
}
int DynFit(const double pars[], int npar, double xtol, double ftol_abs, int MaxIter, int verbose, Estimator_t estimator, Tracer_t * Sample)
{/* MaxLike/MinDist with the given estimator
  output: sample->lnL, the function value at best-fit
		  sample->halo.pars, best-fitting parameters
		  flag_success, return value, whether the search succeeded.
  input: pars, initial guess params
		 npar, number of pars
		 xtol, x-err to consider as convergence
		 ftol_abs, absolute tolerance in function value to consider convergence
		 MaxIter, maximum number of iterations
		 verbose, >0 will switch on verbosed printing
		 estimator, estimator to use
		 sample, tracer sample.
		 */
  const gsl_multimin_fminimizer_type *T =gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *s0, *x0;
  gsl_multimin_function like_func;

  int iter = 0;
  int status, flag_success=0;
  double size, f0,f1;

  /* Starting point */
  x0 = gsl_vector_alloc (npar);
  int i;
  for(i=0;i<npar;i++)  gsl_vector_set (x0, i, pars[i]);

  /* Set initial step sizes to 1 */
  s0 = gsl_vector_alloc (npar);
  gsl_vector_set_all (s0, 1.);
  /* Initialize method and iterate */
  struct like_args args;
  args.estimator=estimator;
  args.Sample=Sample;
  like_func.n = npar;
  like_func.f = like_rfunc;
  like_func.params = &args;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &like_func, x0, s0);

  f1=like_rfunc(x0, &args);
  do {
      iter++;
      int status_iter = gsl_multimin_fminimizer_iterate(s);
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, xtol);
	  f0=f1;
	  f1=s->fval;
	  if (verbose>0)
	  { printf("%5d ", iter);
		for(i=0;i<npar;i++) printf("%g ", gsl_vector_get(s->x, i));
		printf("f()= %g size=%g\n", s->fval, size);
	  }
      if (status == GSL_SUCCESS && fabs(f1-f0)<ftol_abs )
        {
          printf ("Optimization terminated successfully.\n");
		  printf ("\t Current function value: %g\n", f1);
		  printf ("\t Iterations: %d\n", iter);
		  printf ("\t x abs err: %g\n", size);
		  printf ("\t [");
		  for(i=0;i<npar;i++)
			printf("%g ", gsl_vector_get(s->x, i));
		  printf ("\t ]");
		  flag_success=1;
        }
      else if (status_iter)
	  {
		printf("Optimization terminated with error: %d\n", status_iter);
		printf("Warning: %s\n", gsl_strerror (status_iter));
		break;
	  }
    }
  while (status == GSL_CONTINUE && iter < MaxIter);
  
  gsl_vector_free(x0);
  gsl_vector_free(s0);
  gsl_multimin_fminimizer_free (s);

  return flag_success;
}