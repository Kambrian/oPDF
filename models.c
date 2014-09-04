#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
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
double MODEL_TOL_BIN=1e-6, MODEL_TOL_BIN_ABS=1e-6, MODEL_TOL_REL=1e-5; 
//tol_rel=1e-3 is good enough for a contour scan; tol_rel=1e-4 is probably good enough for fmin_gsl() scan; tol_rel=1e-5 should be enough for everything.
//1e-3 should be sufficient, good enough to constrain mass to 1% accuracy with 100000 particles
//accuracy in theta and phase-TS are approximately MODEL_TOL_REL, independent of nP.
//but it's still not accurate enough for minuit to work with the hessian; better use fmin()
double HaloM0,HaloC0,HaloRhos0,HaloRs0,HaloZ0=0.;
int HaloProfID; //profile for interpolation
struct NFWParZ Halo;

struct SplineData
{
	int FlagUseSpline; //spline ready to be used. override the default potential calculation with spline interpolation
	gsl_interp_accel *acc;
	gsl_spline *spline ;
};
static struct SplineData PotSpline;
#pragma omp threadprivate(PotSpline)  //make sure each thread has its own cache
void init_potential_spline()
{
#define LEN_PROF 100
	double PotentialProf[][2][LEN_PROF]={
	  //A4
	  {{1.000000,1.064786,1.133769,1.207222,1.285433,1.368711,1.457384,1.551802,1.652337,1.759385,1.873369,1.994737,2.123968,2.261571,2.408089,2.564099,2.730217,2.907097,3.095436,3.295977,3.509510,3.736877,3.978974,4.236756,4.511238,4.803503,5.114703,5.446064,5.798892,6.174579,6.574605,7.000548,7.454085,7.937005,8.451212,8.998732,9.581724,10.202485,10.863463,11.567263,12.316659,13.114606,13.964249,14.868936,15.832235,16.857942,17.950100,19.113015,20.351270,21.669747,23.073643,24.568491,26.160185,27.854998,29.659611,31.581138,33.627153,35.805721,38.125430,40.595423,43.225437,46.025840,49.007669,52.182678,55.563384,59.163112,62.996052,67.077313,71.422983,76.050190,80.977176,86.223362,91.809427,97.757390,104.090698,110.834316,118.014826,125.660531,133.801572,142.470038,151.700098,161.528137,171.992896,183.135624,195.000244,207.633526,221.085267,235.408492,250.659661,266.898892,284.190198,302.601738,322.206088,343.080524,365.307331,388.974124,414.174193,441.006873,469.577934,500.000000},
	   {218599.539280,218180.855754,217725.429907,217232.192450,216697.537934,216120.716636,215504.340092,214845.000052,214143.042162,213392.693130,212592.951802,211743.032719,210843.471200,209892.037493,208889.560361,207836.104195,206730.081598,205570.840981,204358.391261,203092.946148,201771.227850,200391.672379,198956.574566,197467.823318,195925.984419,194329.100194,192678.451712,190973.428540,189214.425320,187402.715707,185536.534520,183617.295425,181648.085006,179628.275861,177557.937116,175434.126749,173256.908126,171028.196471,168747.973186,166417.164144,164038.090415,161612.679583,159144.165121,156634.643783,154084.755223,151495.932356,148872.818278,146218.836273,143538.649690,140836.859133,138116.288583,135382.806190,132641.144104,129892.828699,127142.712947,124394.558130,121650.169985,118914.842288,116189.370550,113473.635357,110767.296501,108069.003813,105376.815041,102689.301859,100006.909769,97329.544729,94660.935854,92003.470087,89358.758336,86729.791587,84120.219175,81533.130367,78970.712628,76434.064104,73924.603001,71443.335549,68989.406399,66560.206274,64154.743605,61777.961303,59428.191533,57095.154323,54791.916678,52529.117105,50310.985967,48134.577501,46004.280136,43925.864996,41897.462956,39917.960236,37986.735057,36108.395153,34281.850147,32494.075887,30739.450233,29032.798674,27384.318014,25798.139225,24275.334204,22814.149842}
	  },
	  //B4
	  {{1.000000,1.064786,1.133769,1.207222,1.285433,1.368711,1.457384,1.551802,1.652337,1.759385,1.873369,1.994737,2.123968,2.261571,2.408089,2.564099,2.730217,2.907097,3.095436,3.295977,3.509510,3.736877,3.978974,4.236756,4.511238,4.803503,5.114703,5.446064,5.798892,6.174579,6.574605,7.000548,7.454085,7.937005,8.451212,8.998732,9.581724,10.202485,10.863463,11.567263,12.316659,13.114606,13.964249,14.868936,15.832235,16.857942,17.950100,19.113015,20.351270,21.669747,23.073643,24.568491,26.160185,27.854998,29.659611,31.581138,33.627153,35.805721,38.125430,40.595423,43.225437,46.025840,49.007669,52.182678,55.563384,59.163112,62.996052,67.077313,71.422983,76.050190,80.977176,86.223362,91.809427,97.757390,104.090698,110.834316,118.014826,125.660531,133.801572,142.470038,151.700098,161.528137,171.992896,183.135624,195.000244,207.633526,221.085267,235.408492,250.659661,266.898892,284.190198,302.601738,322.206088,343.080524,365.307331,388.974124,414.174193,441.006873,469.577934,500.000000},
	  {109408.840397,109219.265622,109014.953500,108797.046979,108565.648985,108320.268777,108061.349933,107786.842700,107496.639083,107190.273870,106866.830407,106525.719252,106166.577986,105789.315546,105393.545187,104978.910840,104543.439529,104086.509438,103608.436545,103108.857409,102587.835578,102043.970664,101475.265113,100881.969152,100263.273850,99618.967598,98948.030134,98249.333407,97522.023510,96766.538961,95982.123225,95168.044359,94324.017586,93450.260366,92546.015742,91610.142025,90640.923343,89638.414854,88602.810760,87534.087400,86432.618365,85298.431618,84131.050188,82930.851930,81699.057792,80434.460980,79137.847346,77808.462484,76448.028028,75057.418789,73637.862837,72190.858332,70718.048295,69220.655767,67701.771628,66163.140922,64607.485728,63038.715108,61458.438211,59868.885894,58274.273942,56677.353583,55082.151129,53490.955674,51906.741109,50332.571001,48769.315499,47218.915406,45682.518536,44161.074854,42655.165721,41163.355571,39684.252264,38221.793510,36773.532181,35342.491656,33934.766235,32550.426550,31190.062523,29857.506924,28552.385128,27271.166576,26017.185079,24793.231813,23601.002837,22447.689839,21333.133042,20252.987743,19205.196967,18194.444198,17219.782197,16284.065766,15387.277564,14527.046761,13702.833392,12914.898433,12162.582673,11444.065715,10759.436841,10108.487977}
	  }
	};
	if(HaloProfID<0)
	{
	  fprintf(stderr, "Error: HaloProfID=%d, no profile data\n", HaloProfID);
	  exit(1);
	}
	#pragma omp parallel 
	{
		PotSpline.acc= gsl_interp_accel_alloc ();
		PotSpline.spline= gsl_spline_alloc (gsl_interp_cspline, LEN_PROF);
		gsl_spline_init(PotSpline.spline, PotentialProf[HaloProfID][0], PotentialProf[HaloProfID][1], LEN_PROF);	
		PotSpline.FlagUseSpline=1;
	}
}
void free_potential_spline()
{
	#pragma omp parallel
	{
	 if(PotSpline.FlagUseSpline)
	 {
	  gsl_spline_free (PotSpline.spline);
	  gsl_interp_accel_free(PotSpline.acc);
	  PotSpline.FlagUseSpline=0;
	 }
	}
}
double eval_potential_spline(double r)
{
  return -gsl_spline_eval(PotSpline.spline, r, PotSpline.acc); //the input are |pot|, restore the sign.
}

void define_halo(const double pars[])
{ 
#if FIT_PAR_TYPE==PAR_TYPE_M_C
  decode_NFWprof(HaloZ0,pars[0]*HaloM0,pars[1]*HaloC0,VIR_C200,&Halo);  
#elif  FIT_PAR_TYPE==PAR_TYPE_RHOS_RS
  decode_NFWprof2(HaloZ0,pars[0]*HaloRhos0,pars[1]*HaloRs0,VIR_C200,&Halo);
#elif FIT_PAR_TYPE==PAR_TYPE_POTS_RS  
  decode_NFWprof2(HaloZ0,pars[0]/pars[1]/pars[1]*HaloRhos0,pars[1]*HaloRs0,VIR_C200,&Halo);
// #elif FIT_PAR_TYPE==PAR_TYPE_INTERPOLATE
//   init_potential_spline();
#endif
}

double NFW_mass(double r)
{ 
  double x=r/Halo.Rs;
  return Halo.Ms*(log(1+x)-x/(1+x)); 
}
double NFW_like(double pars[], Tracer_t *T)
{
  if(pars[0]<=0||pars[1]<=0||isnan(pars[0])||isnan(pars[1])) return -INFINITY;
  define_halo(pars);
  double lnL=log(Halo.Rhos)*T->nP-(NFW_mass(T->rmax)-NFW_mass(T->rmin))/T->mP; //the normalizations
  int i;
  #pragma omp parallel for reduction(+: lnL)
  for(i=0;i<T->nP;i++)
  {
	double r=T->P[i].r/Halo.Rs;
	lnL+=-log(r)-2.*log(1.+r);
  }
  return lnL; //loglikelihood
}

double halo_pot(double r)
{
  if(PotSpline.FlagUseSpline)  return eval_potential_spline(r); //use spline if inited.
  
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
	  if (dx/( P->r-lbin[0])<MODEL_TOL_BIN||dx<MODEL_TOL_BIN_ABS) break;
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
	  if ( is_forbidden ( xmid, P->E, P->L2 ) )
		rbin[1]=xmid;
	  else
		rbin[0]=xmid;
	  if ( dx/(rbin[1]-P->r ) <MODEL_TOL_BIN||dx<MODEL_TOL_BIN_ABS ) break;
	}
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
  //   gsl_integration_qags (&F, xlim[0],xlim[1], 0, MODEL_TOL_REL, MODEL_MAX_INTVAL, //3, 
  // 		       GSL_workspace, &result, &error);
  size_t neval;
  gsl_integration_cquad (&F, P->rlim[0], P->rlim[1], 0, MODEL_TOL_REL, //3, 
			 GSL_workspaceC, &(P->T), &error, &neval);
  if(P->T<=0||isnan(P->T)||P->T==INFINITY){ 
	fprintf(stderr,"Warning: T=%g, reset to 1. [%g,%g]\n", P->T, P->rlim[0], P->rlim[1]);P->T=1.;}
  if(IS_PHASE_ESTIMATOR(estimator))
  {
  double t;
  gsl_integration_cquad (&F, P->rlim[0],P->r, 0, MODEL_TOL_REL, //3, 
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
	  if(Sample->FlagUseWeight)
		c+=Sample->P[i].theta*Sample->P[i].w;
	  else
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
	  if(Sample->FlagUseWeight)
		c+=Sample->P[i].theta*Sample->P[i].w;
	  else
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
  printf("Old Test\n");
  #pragma omp parallel
  {
   #pragma omp for 
    for(i=0;i<Sample->nP;i++)
      theta[i]=Sample->P[i].theta;
    #pragma omp single
      {
      qsort(theta, Sample->nP, sizeof(double), cmpDouble);
      lnL=-log((1-theta[0])*theta[Sample->nP-1])-theta[0]+theta[Sample->nP-1]-1; //head and tail partitions
      printf("head: %g\n", lnL);
      }
    #pragma omp for reduction(+:lnL)
    for(i=1;i<Sample->nP;i++)
	  lnL+=log(theta[i]/theta[i-1])*i*i/Sample->nP/Sample->nP
	  -log((1-theta[i])/(1-theta[i-1]))*(1.-(double)i/Sample->nP)*(1.-(double)i/Sample->nP)
	  +theta[i-1]-theta[i];
//		lnL+=log(P[i].theta/P[i-1].theta)*P[i-1].cw*P[i-1].cw/Sample->nP/Sample->nP
// 	  -log((1-P[i].theta)/(1-P[i-1].theta))*(1.-P[i-1].cw/Sample->nP)*(1.-P[i-1].cw/Sample->nP)
// 	  +P[i-1].theta-P[i].theta;
  }
  free(theta);
  lnL*=-Sample->nP; //differ by -1, to make it comparable to a loglike
  return lnL;
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
static double ADgpar[2][3]={{0.569,   -0.570,    0.511}, { 0.431,    0.227,    0.569}}; //w, mu, sigma
double AndersonDarlingLike(Tracer_t *Sample, int estimator)
{/* adopting a PDF fit of ln(AD),
    return the log(Prob), the real likelihood (posterior probability)
    */
    double lnAD=log(-AndersonDarlingTest(Sample));
	switch(estimator)
	{
	  case RADIAL_PHASE_AD_GEV:
		return log(GeneralizedExtremeValuePDF(lnAD, LnAD_GEV_MU, LnAD_GEV_SIGMA, LnAD_GEV_K)); //loglikelihood
	  case RADIAL_PHASE_AD_BINORMAL:
		return log(ADgpar[0][0]*NormPDF(lnAD, ADgpar[0][1], ADgpar[0][2])
				  +ADgpar[1][0]*NormPDF(lnAD, ADgpar[1][1], ADgpar[1][2]));
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

void predict_radial_count (double RadialCountPred[], int nbin, int FlagRLogBin, Tracer_t *Sample)
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
		F.params = &Fpar;
		gsl_integration_cquad (&F, MAX(rbin[0],Sample->P[j].rlim[0]), MIN(rbin[1],Sample->P[j].rlim[1]), 0, MODEL_TOL_REL, //3, 
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

void like_init(const double pars[], int estimator, Tracer_t *Sample)
{
  int i;
  if(pars[0]<=0||pars[1]<=0||isnan(pars[0])||isnan(pars[1])) return;
  define_halo(pars);
  
  #pragma omp parallel for
  for(i=0;i<Sample->nP;i++)
    solve_radial_orbit(Sample->P+i,Sample->rmin,Sample->rmax,estimator);
}
double like_eval(const double pars[], int estimator,Tracer_t *Sample)
{
  double lnL;
  
  if(pars[0]<=0||pars[1]<=0||isnan(pars[0])||isnan(pars[1])) return -INFINITY;
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
	case (RADIAL_PHASE_ROULETTE+100):
      lnL=AndersonDarlingTest_Old(Sample);
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
  Sample->lnL=lnL;
  return lnL;
}

double likelihood(const double pars[], int estimator, Tracer_t *Sample)
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
double wenting_like(const double pars[], Tracer_t *Sample)
{//pars in units of real parameters.
  int i;
  if(pars[0]<=0||pars[1]<=0||isnan(pars[0])||isnan(pars[1])) return -INFINITY;
  define_halo(pars);
  
  double wpar[NUM_PAR_MAX],lnL=0.;
  wpar[0]=Halo.Rhos*1e10; //mass unit Msun
  wpar[1]=Halo.Rs;
  wpar[2]=pars[2]*0.715;
  wpar[3]=pars[3]*69.014;
  wpar[4]=pars[4]*2.301;
  wpar[5]=pars[5]*7.467;
  #pragma omp parallel for reduction(+:lnL) 
  for(i=0;i<Sample->nP;i++)
	lnL+=dataprob_model(Sample->P[i].r,Sample->P[i].vr,sqrt(Sample->P[i].L2)/Sample->P[i].r, wpar);
//   printf("%g,%g,%g,%g,%g,%g: %g\n", pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],lnL);
  Sample->lnL=lnL;
  return lnL; //loglikelihood
}

double like_to_chi2(double lnL, int estimator)
{//convert likelihood() values to a chi-square measure
  double lnAD;
  switch(estimator)
  {
	case RADIAL_PHASE_LMEANRAW:
	  return lnL*lnL;
	case RADIAL_PHASE_LMOMENT:
	  return -lnL;
	case RADIAL_PHASE_ROULETTE:
	  return -lnL; //directly use AD as a chi-2
// 	  lnAD=log(-lnL);
// 	  return  -2.*log(ADgpar[0][0]*NormPDF(lnAD, ADgpar[0][1], ADgpar[0][2])
// 				  +ADgpar[1][0]*NormPDF(lnAD, ADgpar[1][1], ADgpar[1][2]));
// 	  lnL=(log(-lnL)-LnADMean)/LnADSig;
// 	  return lnL*lnL;
	case RADIAL_PHASE_AD_GEV:
	case RADIAL_PHASE_AD_BINORMAL:
	case RADIAL_PHASE_AD_NORMAL:  
	case RADIAL_BIN_ESTIMATOR:
	  return -2.*lnL; //simply 2*(the negative loglike), to make it comparable to chisquare; Note these could differ from the chi-2 by a constant.
	default:
	  fprintf(stderr, "Error: like_to_chi2() not implemented for estimator=%d yet\n", estimator);
	  exit(1);
  }
}

void freeze_energy(const double pars[], Tracer_t *Sample)
{//fix the energy parameter according to initial potential
  int i;
  if(pars[0]<=0||pars[1]<=0||isnan(pars[0])||isnan(pars[1])) return;
  define_halo(pars);
  #pragma omp parallel for
  for(i=0;i<Sample->nP;i++)
      Sample->P[i].E=-(Sample->P[i].K+halo_pot(Sample->P[i].r));//differ from previous version
}

double freeze_and_like(const double pars[], int estimator, Tracer_t *Sample)
{
//   printf("%g,%g\n", pars[0], pars[1]);
  freeze_energy(pars, Sample);
  return likelihood(pars, estimator, Sample);
}

double jointLE_FChi2(const double pars[], int estimator, int nbinL, int nbinE, Tracer_t *Sample)
{//automatically update views if needed; then freeze and like
  int i;
  double chi2;
  if(Sample->nView!=nbinL||Sample->ViewType!='L') //already allocated
	create_tracer_views(Sample, nbinL, 'L');
  for(i=0,chi2=0;i<nbinL;i++)
	chi2+=jointE_FChi2(pars, estimator, nbinE, Sample->Views+i);
  return chi2;
}
double jointE_FChi2(const double pars[], int estimator, int nbin, Tracer_t *Sample)
{//this does freeze_and_like
  int i;
  double lnL, chi2;
  freeze_energy(pars, Sample);
  create_tracer_views(Sample, nbin, 'E');
  for(i=0,chi2=0;i<nbin;i++)
  {
	if(estimator==RADIAL_BIN_ESTIMATOR)  count_tracer_radial(Sample->Views+i, NumRadialCountBin, 1);
	lnL=likelihood(pars, estimator, Sample->Views+i);
	chi2+=like_to_chi2(lnL, estimator);
	if(estimator==RADIAL_BIN_ESTIMATOR)  free_tracer_rcounts(Sample->Views+i);
  }
  free_tracer_views(Sample);
  Sample->lnL=chi2;
  return chi2;
}
void create_nested_views(const double pars[], int nbin[], char ViewTypes[], Tracer_t *Sample)//iteratively: freeze, fit; then freeze, then fit.
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
double nested_views_Chi2(const double pars[], int estimator, Tracer_t *Sample)
{//pure like, without freezing energy; freeze it yourself before calling this!
  double lnL;
  if(!Sample->nView)
  {
	lnL=likelihood(pars, estimator, Sample);
	return like_to_chi2(lnL,estimator);
  }
  
  int i;
  for(i=0,lnL=0.;i<Sample->nView;i++)
	lnL+=nested_views_Chi2(pars, estimator, Sample->Views+i);
  Sample->lnL=lnL;
  return lnL;
}
double nested_views_FChi2(const double pars[], int estimator, Tracer_t *Sample)
{//freeze and like
  freeze_energy(pars, Sample);
  return nested_views_Chi2(pars, estimator, Sample);
}