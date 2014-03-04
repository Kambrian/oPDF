#include "cosmology.h"

#define NUM_PAR_MAX 10

extern struct NFWParZ Halo;

extern void alloc_integration_space();
extern void free_integration_space();
extern double halo_pot(double r);
extern void solve_radial_limits(int pid);
extern double vr_inv_part(double r, int pid);
extern void solve_radial_orbit(int pid);
extern double likelihood(double pars[]);
