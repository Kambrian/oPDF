#include "globals.h"
struct global Globals;

void default_global_pars()
{
  Globals.tol.bin=1e-6;
  Globals.tol.bin_abs=1e-6;
  Globals.tol.rel=1e-5;
  //tol_rel=1e-3 is good enough for a contour scan; tol_rel=1e-4 is probably good enough for fmin_gsl() scan; tol_rel=1e-5 should be enough for everything.
//1e-3 should be sufficient, good enough to constrain mass to 1% accuracy with 100000 particles
//accuracy in theta and phase-TS are approximately MODEL_TOL_REL, independent of nP.
//but it's still not accurate enough for minuit to work with the hessian; better use fmin()
  Globals.cosmology.OmegaM0=0.3;
  Globals.cosmology.OmegaL0=0.7;
  set_units(1e10, 1., 1.);//default units set to 1e10Msun/h, kpc/h, km/s
}

void set_units(double MassInMsunh, double LengthInKpch, double VelInKms)
{//all functions below in internal units, determined by set_units()
  Globals.units.MassInMsunh=MassInMsunh;
  Globals.units.LengthInKpch=LengthInKpch;
  Globals.units.VelInKms=VelInKms;
  Globals.units.Const.G=4.30071e-6*MassInMsunh/VelInKms/VelInKms/LengthInKpch;
  Globals.units.Const.H0=0.1*LengthInKpch/VelInKms;
}