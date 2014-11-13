#ifndef MODELS_HEADER_INCLUDED

#define HALF_ORBIT_PERIOD 1
#define FULL_ORBIT_PERIOD 2
#define PHASE_PERIOD HALF_ORBIT_PERIOD
// #define SWAP_T
typedef enum
{
  EID_RBinLike=0,
  EID_PhaseAD,
  EID_PhaseMean,
  EID_PhaseMeanRaw
} Estimator_t;

extern void alloc_integration_space();
extern void free_integration_space();

extern void solve_radial_limits ( Particle_t *P, double rmin, double rmax, Halo_t *halo);
extern double vr_inv_part(double r, double E, double L2, Halo_t *halo);
extern void solve_radial_orbit(Particle_t *P, double rmin, double rmax, Halo_t *halo, int FlagSetPhase);

extern void predict_radial_count (double RadialCountPred[], int nbin, int FlagRLogBin, Tracer_t *Sample);

extern double like_eval(Estimator_t estimator, Tracer_t *Sample);
extern double nested_views_like(Estimator_t estimator,  Tracer_t *Sample, int nbin_r, int FlagRLogBin);
extern int DynFit(const double pars[], int npar, double xtol, double ftol_abs, int MaxIter, int verbose, Estimator_t estimator, Tracer_t * Sample);
#define MODELS_HEADER_INCLUDED
#endif