#include "cosmology.h"

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

#ifndef NBIN_R
#define NBIN_R 30
#endif

//#define ESTIMATOR 10  //this has been moved to makefile
//#define RETURN_RAWMEAN 
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

#define MODEL_TOL_BIN 1e-6 //bin size relative error
#define MODEL_TOL_REL 1e-8 //relative tolerance for PERIOD integral
#define MODEL_MAX_INTVAL 1000

extern double M0,C0,Rhos0,Rs0;
extern struct NFWParZ Halo;

extern void alloc_integration_space();
extern void free_integration_space();
extern double halo_pot(double r);
extern void solve_radial_limits(int pid);
extern double vr_inv_part(double r, int pid);
extern void solve_radial_orbit(int pid, int estimator);
extern void fill_radial_bin();
extern double likelihood(double pars[], int estimator);

extern void define_halo(double pars[]);
extern void init();
extern void freeze_energy(double pars[]);
extern double freeze_and_like(double pars[], int estimator);
extern void select_particles(int subsample_id);