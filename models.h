#include "cosmology.h"

#define NUM_PAR_MAX 10

extern struct NFWParZ Halo;

extern void alloc_integration_space();
extern void free_integration_space();
extern double halo_pot(double r);
extern double vr_inv(double r, void *params);
extern double Period(double r, double K, double L2, double rmin, double rmax);
extern double likelihood(double pars[]);
