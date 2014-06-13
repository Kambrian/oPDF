#include "cosmology.h"

#define PAR_TYPE_M_C 0
#define PAR_TYPE_RHOS_RS 1
#define FIT_PAR_TYPE 0

#define HALF_ORBIT_PERIOD 1
#define FULL_ORBIT_PERIOD 2

#define PHASE_PERIOD 1
// #define SWAP_T

// typedef enum
// {
//   EID_Radial=0,
//   EID_RadialBin,
//   EID_MixedRadial,
//   EID_Entropy,
//   EID_Iterative,
//   EID_Phase=100,
//   EID_PhaseAD,
//   EID_PhaseADProb,
//   EID_PhaseMean,
//   EID_PhaseMeanRaw,
//   EID_PhaseBin,
//   EID_PhasePartition,
//   EID_PhaseProcess,
//   EID_PhaseKS,
//   EID_PhaseKuiper,
//   EID_PhaseResultant,
//   EID_PhaseCosMean,
// } Estimator_t;
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

//#define ESTIMATOR 10  //this has been moved to makefile
// #define RETURN_RAWMEAN 
// #define RETURN_PROB //for KS and Kuiper

#define IS_PHASE_ESTIMATOR(x) ((x)>=RADIAL_PHASE_BINNED)

#ifndef PHASE_PERIOD
  #if ESTIMATOR==RADIAL_PHASE_KUIPER||ESTIMATOR==RADIAL_PHASE_CMOMENT||ESTIMATOR==RADIAL_PHASE_COSMEAN //these two seems to be biased a bit when used with the half-period-phase
    #ifndef SWAP_T 
      #define PHASE_PERIOD FULL_ORBIT_PERIOD
    #else //alternative period def
      #define PHASE_PERIOD HALF_ORBIT_PERIOD
    #endif
  #else //below can be biased if used with the full-period-phase
    #ifndef SWAP_T
      #define PHASE_PERIOD HALF_ORBIT_PERIOD 
    #else
      #define PHASE_PERIOD FULL_ORBIT_PERIOD 
    #endif
  #endif
#endif


#define NUM_PAR_MAX 10
extern double MODEL_TOL_BIN, MODEL_TOL_BIN_ABS, MODEL_TOL_REL;
extern double HaloM0,HaloC0,HaloRhos0,HaloRs0,HaloZ0; //to define the reference point (units) of the scan
extern struct NFWParZ Halo;

extern void alloc_integration_space();
extern void free_integration_space();

extern void define_halo(double pars[]);
extern double halo_pot(double r);
extern void solve_radial_limits ( Particle_t *P, double rmin, double rmax);
extern double vr_inv_part(double r, double E, double L2);
extern void solve_radial_orbit(Particle_t *P, double rmin, double rmax, int estimator);

extern void predict_radial_count (double RadialCountPred[], int nbin, Tracer_t *Sample);

extern void like_init(double pars[], int estimator, Tracer_t *Sample);
extern double like_eval(double pars[], int estimator,Tracer_t *Sample);
extern double likelihood(double pars[], int estimator, Tracer_t *Sample);
extern double like_to_chi2(double lnL, int estimator);

extern void freeze_energy(double pars[], Tracer_t *Sample);
extern double freeze_and_like(double pars[], int estimator, Tracer_t *Sample);

extern double jointLE_like(double pars[], int estimator, int nbinL, int nbinE, Tracer_t *Sample);
extern double jointE_like(double pars[], int estimator, int nbin, Tracer_t *Sample);