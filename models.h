#include "cosmology.h"

#define MIXED_RADIAL_ESTIMATOR 0
#define ITERATIVE_ESTIMATOR 1
#define ENTROPY_ESTIMATOR 2

#define ESTIMATOR MIXED_RADIAL_ESTIMATOR

#define NUM_PAR_MAX 10

#define MODEL_TOL_BIN 1e-3 //bin size relative error
#define MODEL_TOL_REL 1e-6 //relative tolerance for PERIOD integral
#define MODEL_MAX_INTVAL 1000

extern struct NFWParZ Halo;

extern void alloc_integration_space();
extern void free_integration_space();
extern double halo_pot(double r);
extern void solve_radial_limits(int pid);
extern double vr_inv_part(double r, int pid);
extern void solve_radial_orbit(int pid);
extern double likelihood(double pars[]);

extern void init();
extern void freeze_energy(double pars[]);
extern double freeze_and_like(double pars[]);