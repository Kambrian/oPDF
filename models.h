#ifndef MODELS_HEADER_INCLUDED
#define MODELS_HEADER_INCLUDED

#include "cosmology.h"
#include "wenting.h"

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

#define LnADMean (-0.22)
#define LnADSig 0.66
#define LnAD_GEV_K (-0.2128)
#define LnAD_GEV_SIGMA 0.6335
#define LnAD_GEV_MU (-0.4776)

#define IS_PHASE_ESTIMATOR(x) ((x)>=EID_PhaseAD)

extern double MODEL_TOL_BIN, MODEL_TOL_BIN_ABS, MODEL_TOL_REL;
extern struct NFWParZ Halo;

//potential profile using interpolation
extern void init_potential_spline();
extern void free_potential_spline();
extern double eval_potential_spline(double r);
extern double eval_density_spline(double r);

extern void alloc_integration_space();
extern void free_integration_space();

extern void decode_TemplateProf(double z, double M, double c, int virtype, struct NFWParZ *halo);
extern void decode_TemplateProf2(double z, double Pots, double Rs, int virtype, struct NFWParZ *halo);

extern double NFW_like(double pars[], Tracer_t *T);
extern void define_halo(const double pars[]);
extern double halo_pot(double r);
extern double halo_mass(double r);
extern void solve_radial_limits ( Particle_t *P, double rmin, double rmax);
extern double vr_inv_part(double r, double E, double L2);
extern void solve_radial_orbit(Particle_t *P, double rmin, double rmax, int estimator);

extern void predict_radial_count (double RadialCountPred[], int nbin, int FlagRLogBin, Tracer_t *Sample);

/* --types of likelihood:
likelihood: elemetary likelihood evaluation =(init+eval)
freeze_and_like: freeze energy and calc like =(freeze_energy+likelihood)
joint_FChi2: automatically create/manage views, freeze energy, and likelihood
views_Chi2: likelihood on existing views
views_FChi2: freeze_and_like on existing views
**_Chi2 functions returns -2*log(like).
*/

extern void like_init(const double pars[], int estimator, Tracer_t *Sample); //define halo and prepare orbit
extern double like_eval(const double pars[], int estimator,Tracer_t *Sample); 
extern double likelihood(const double pars[], int estimator, Tracer_t *Sample); //init and eval
extern double wenting_like(const double pars[], Tracer_t *Sample);
extern double like_to_chi2(double lnL, int estimator);

extern void freeze_energy(const double pars[], Tracer_t *Sample);
extern double freeze_and_like(const double pars[], int estimator, Tracer_t *Sample);

extern double jointLE_FChi2(const double pars[], int estimator, int nbinL, int nbinE, Tracer_t *Sample);
extern double jointE_FChi2(const double pars[], int estimator, int nbin, Tracer_t *Sample);

extern void create_nested_views(const double pars[], int nbin[], char ViewTypes[], Tracer_t *Sample);
extern double nested_views_Chi2(const double pars[], int estimator, Tracer_t *Sample);
extern double nested_views_FChi2(const double pars[], int estimator, Tracer_t *Sample);

#endif