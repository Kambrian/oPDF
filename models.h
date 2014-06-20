#include "cosmology.h"

#define NR_INTEGRAL 1
#define GSL_INTEGRAL 2
extern int NumericalIntegralRoutines, RLimValConv, RLimInnerVal;
// #define NR_INTEGRAL //use NumericalRecipes integal instead of GSL.

#define PAR_TYPE_M_C 0
#define PAR_TYPE_RHOS_RS 1
#define FIT_PAR_TYPE 0

#define HALF_ORBIT_PERIOD 1
#define FULL_ORBIT_PERIOD 2

#define PHASE_PERIOD 1
// #define SWAP_T

#define ENTROPY_ESTIMATOR 1 //conditional like; junk
#define ITERATIVE_ESTIMATOR 2 //iterative entropy; still junk
#define MIXED_RADIAL_ESTIMATOR 3
#define RADIAL_BIN_ESTIMATOR 4
#define RADIAL_PHASE_BINNED 5
#define RADIAL_PHASE_PARTITION 6//multinomial distribution
#define RADIAL_PHASE_PROCESS 7    //Poisson process time delay; junk
#define RADIAL_PHASE_ROULETTE 8
#define RADIAL_PHASE_CMOMENT 9  //circular moment
#define RADIAL_PHASE_LMOMENT 10 //linear moment, under roulette phase definition
#define RADIAL_PHASE_KS 11 //ks-test
#define RADIAL_PHASE_KUIPER 12 //kuiper's test
#define RADIAL_PHASE_COSMEAN 13 //<cos(theta)>
#define RADIAL_PHASE_LMEANRAW 14 //mean phase, std normal
#define RADIAL_PHASE_AD_GEV 15 //ADtest likelihood, from GEV fit
#define RADIAL_PHASE_AD_BINORMAL 16 //AD like, from binormal fit
#define RADIAL_PHASE_AD_NORMAL 17 //AD like, from normal fit
#define LnADMean (-0.22)
#define LnADSig 0.66
#define LnAD_GEV_K (-0.2128)
#define LnAD_GEV_SIGMA 0.6335
#define LnAD_GEV_MU (-0.4776)

// #define RETURN_PROB //for KS and Kuiper

#define IS_PHASE_ESTIMATOR(x) ((x)>=RADIAL_PHASE_BINNED)


#define NUM_PAR_MAX 10
extern double MODEL_TOL_BIN, MODEL_TOL_BIN_ABS, MODEL_TOL_REL;
extern double HaloM0,HaloC0,HaloRhos0,HaloRs0,HaloZ0; //to define the reference point (units) of the scan
extern struct NFWParZ Halo;

extern void choose_integral_routines(int type, int rlimvalconv, int rliminnerval);
extern void alloc_integration_space();
extern void free_integration_space();

extern void define_halo(double pars[]);
extern double halo_pot(double r);
extern void solve_radial_limits ( Particle_t *P, double rmin, double rmax);
extern double vr_inv_part(double r, double E, double L2);
extern void solve_radial_orbit(Particle_t *P, double rmin, double rmax, int estimator);

extern void predict_radial_count (double RadialCountPred[], int nbin, Tracer_t *Sample);

/* --types of likelihood:
likelihood: elemetary likelihood evaluation =(init+eval)
freeze_and_like: freeze energy and calc like =(freeze_energy+likelihood)
joint_Flike: automatically create/manage views, freeze energy, and likelihood
views_like: likelihood on existing views
views_Flike: freeze_and_like on existing views
*/

extern void like_init(double pars[], int estimator, Tracer_t *Sample); //define halo and prepare orbit
extern double like_eval(double pars[], int estimator,Tracer_t *Sample); 
extern double likelihood(double pars[], int estimator, Tracer_t *Sample); //init and eval
extern double like_to_chi2(double lnL, int estimator);

extern void freeze_energy(double pars[], Tracer_t *Sample);
extern double freeze_and_like(double pars[], int estimator, Tracer_t *Sample);

extern double jointLE_Flike(double pars[], int estimator, int nbinL, int nbinE, Tracer_t *Sample);
extern double jointE_Flike(double pars[], int estimator, int nbin, Tracer_t *Sample);

extern void create_nested_views(double pars[], int nbin[], char ViewTypes[], Tracer_t *Sample);
extern double nested_views_like(double pars[], int estimator, Tracer_t *Sample);
extern double nested_views_Flike(double pars[], int estimator, Tracer_t *Sample);
