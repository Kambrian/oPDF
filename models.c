#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mymath.h"
#include "globals.h"
#include "cosmology.h"
#include "tracer.h"
#include "models.h"
#include "halo.h"
	
#define MODEL_MAX_INTVAL 1000

double NFW_like(double pars[], Tracer_t *T, Halo_t *halo)
{//define the appropriate halo_type first
  if(pars[0]<=0||pars[1]<=0||isnan(pars[0])||isnan(pars[1])) return -INFINITY;
  halo_set_param(pars, halo);
  double lnL=log(halo->Rhos)*T->nP-(halo_mass(T->rmax, halo)-halo_mass(T->rmin, halo))/T->mP; //the normalizations
  int i;
  #pragma omp parallel for reduction(+: lnL)
  for(i=0;i<T->nP;i++)
  {
	double r=T->P[i].r/halo->Rs;
	lnL+=-log(r)-2.*log(1.+r);
  }
  return lnL; //loglikelihood
}

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
  double vr2=2*(-E-L2/2./r/r-halo_pot(r, halo));//what about hubble flow?
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
void solve_radial_orbit(Particle_t *P, double rmin, double rmax, int estimator, Halo_t *halo)//bottleneck in gsl_integration_cquad
{//find peri(apo)-centers and integrate the period 
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
  if(IS_PHASE_ESTIMATOR(estimator))
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

double like_mean_phase_raw(Tracer_t *Sample)
{  //raw mean phase (not squared)
  int i;
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<Sample->nP;i++)
    {
	  if(Sample->FlagUseWeight)
		c+=Sample->P[i].theta*Sample->P[i].w;
	  else
		c+=Sample->P[i].theta;
    }
    c/=Sample->nP;c-=0.5; 
    return c*sqrt(12.*Sample->nP);//standard normal variable
}
double like_mean_phase(Tracer_t *Sample)
{ //mean phase 
  int i;
  double c=0.;
   #pragma omp parallel for reduction(+:c)
    for(i=0;i<Sample->nP;i++)
    {
	  if(Sample->FlagUseWeight)
		c+=Sample->P[i].theta*Sample->P[i].w;
	  else
		c+=Sample->P[i].theta;
    }
    c/=Sample->nP;c-=0.5; 
  //squared mean phase
    return -c*c*12*Sample->nP; //a standard chisquare variable
}
typedef struct
{
  double theta;
  double w;//weight
  double cw;//cumsum of w
} ThetaWeight_t;
static int cmpThetaWeight(const void *p1, const void *p2)
{ //in ascending order
  if(((ThetaWeight_t *)p1)->theta > ((ThetaWeight_t *)p2)->theta ) 
    return 1;
  
  if(((ThetaWeight_t *)p1)->theta < ((ThetaWeight_t *)p2)->theta )
    return -1;
  
  return 0;
}
double AndersonDarlingTest(Tracer_t *Sample)
{//Beloborodov&Levin 2004, apj, 613:224-237; simplified equation as in the note.
  int i;
  double lnL=0.;
  ThetaWeight_t *P;
  P=malloc(sizeof(ThetaWeight_t)*Sample->nP);
  
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<Sample->nP;i++)
	{
      P[i].theta=Sample->P[i].theta;
	  P[i].w=Sample->P[i].w;
	}
    #pragma omp single
    {
      qsort(P, Sample->nP, sizeof(ThetaWeight_t), cmpThetaWeight);
	  if(Sample->FlagUseWeight)
	  {
		double wsum=0.;
		for(i=0;i<Sample->nP;i++)
		{
		  wsum+=P[i].w;
		  P[i].cw=wsum;
		}
	  }
	}
    #pragma omp for reduction(+:lnL)
    for(i=0;i<Sample->nP;i++)
	{
	  if(Sample->FlagUseWeight)
		lnL+=P[i].w*(P[i].w-2.*P[i].cw)*log(P[i].theta)-P[i].w*(P[i].w+2.*(Sample->nP-P[i].cw))*log(1-P[i].theta);
	  else
		lnL+=-(1.+2.*i)*log(P[i].theta)+(1+2.*(i-Sample->nP))*log(1-P[i].theta); //note the i starts from 0 here, so they are actually j=i-1, with i starting from 1
	}
  }
  free(P);
  lnL=lnL/Sample->nP-Sample->nP; 
  return -lnL; //differ by -1, to make it comparable to a loglike
}
void predict_radial_count (double RadialCountPred[], int nbin, int FlagRLogBin,  Tracer_t *Sample)
{
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
			Fpar.halo=halo;
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
		Fpar.halo=halo;
		F.params = &Fpar;
		gsl_integration_cquad (&F, MAX(rbin[0],Sample->P[j].rlim[0]), MIN(rbin[1],Sample->P[j].rlim[1]), 0, Globals.tol.rel, //3, 
			       GSL_workspaceC, &t, &error, &neval);
		if(Sample->FlagUseWeight)
		  p+=t/Sample->P[j].T*Sample->P[j].w;
		else
		  p+=t/Sample->P[j].T;
      }
      lnL+=Sample->RadialCount[i]*log(p);
    }
    
    return lnL;
}

void tracer_set_orbits(int estimator,  Tracer_t *Sample)
{
  int i;

  #pragma omp parallel for
  for(i=0;i<Sample->nP;i++)
    solve_radial_orbit(Sample->P+i,Sample->rmin,Sample->rmax,estimator,Sample->halo);
}

double like_eval(int estimator, Tracer_t *Sample)
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
      lnL=like_mean_phase(Sample);
      break;
	case EID_PhaseMeanRaw:
      lnL=like_mean_phase_raw(Sample);
      break;  
    default:
      DEBUGPRINT("Error: unknown Estimator=%d\n", estimator);
      exit(estimator);
  }
//     printf("%g,%g: %g\n", pars[0],pars[1],lnL);
  Sample->lnL=lnL;
  return lnL;
}

double likelihood(int estimator, Tracer_t *Sample)
{
  tracer_set_orbits(estimator, Sample);
  double lnL=like_eval(estimator, Sample);
  return lnL;
}

double jointLE_FChi2(const double pars[], int estimator, int nbinL, int nbinE,  Tracer_t *Sample)
{//automatically update views if needed; then freeze and like
  int i;
  double chi2;
  if(Sample->nView!=nbinL||Sample->ViewType!='L') //already allocated
	create_tracer_views(Sample, nbinL, 'L');
  for(i=0,chi2=0;i<nbinL;i++)
	chi2+=jointE_FChi2(pars, estimator, nbinE, halo, Sample->Views+i);
  return chi2;
}
double jointE_FChi2(const double pars[], int estimator, int nbin,  Tracer_t *Sample)
{//this does freeze_and_like
  int i;
  double lnL, chi2;
  freeze_energy(halo, Sample);
  create_tracer_views(Sample, nbin, 'E');
  for(i=0,chi2=0;i<nbin;i++)
  {
	if(estimator==EID_RBinLike)  count_tracer_radial(Sample->Views+i, NumRadialCountBin, 1);
	lnL=likelihood(estimator, halo, Sample->Views+i);
	chi2+=like_to_chi2(lnL, estimator);
	if(estimator==EID_RBinLike)  free_tracer_rcounts(Sample->Views+i);
  }
  free_tracer_views(Sample);
  Sample->lnL=chi2;
  return chi2;
}
void create_nested_views(int nbin[], char ViewTypes[],  Tracer_t *Sample)
{//ViewTypes should be a non-empty string!
  //used to create static views (jointE, jointLE dynamically create views by themselves instead)
  //do not mix this with joint_like functions, since they destroy the views!
  int i;
  if(ViewTypes[0]=='E')  freeze_energy(halo, Sample);
  create_tracer_views(Sample, nbin[0], ViewTypes[0]);
  if(ViewTypes[1]=='\0') return; //done
  
  for(i=0;i<nbin[0];i++)
	create_nested_views(nbin+1, ViewTypes+1, halo, Sample->Views+i);//descend
}
double nested_views_Chi2(int estimator,  Tracer_t *Sample)
{//pure like, without freezing energy; freeze it yourself before calling this!
  double lnL;
  if(!Sample->nView)
  {
	lnL=likelihood(estimator, Sample);
	return like_to_chi2(lnL,estimator);
  }
  
  int i;
  for(i=0,lnL=0.;i<Sample->nView;i++)
	lnL+=nested_views_Chi2(estimator, halo, Sample->Views+i);
  Sample->lnL=lnL;
  return lnL;
}
double nested_views_FChi2(int estimator,  Tracer_t *Sample)
{//freeze and like
  freeze_energy(halo, Sample);
  return nested_views_Chi2(estimator, halo, Sample);
}