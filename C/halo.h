#ifndef HALO_HEADER_INCLUDED

// #include "globals.h"

#define NUM_PAR_MAX 10

typedef enum
{
  HT_NFWMC=0,
  HT_NFWPotsRs,
  HT_NFWRhosRs,
  HT_TMPMC,
  HT_TMPPotScaleRScale,
  HT_CoreRhosRs,
  HT_CorePotsRs,
  HT_PointM,
  HT_IsothermalK
} HaloType_t;

typedef struct Halo
{//all quantities are physical
  double pars[NUM_PAR_MAX]; //raw parameters
  double scales[NUM_PAR_MAX]; //parameter scales
  double z;
  double M;
  double c;
  double Rv;
  double Pots;//-4*pi*G*rhos*rs^2, the potential at r=0
  double Rs;
  double Rhos;
  double Ms;//4*pi*rs^3*rhos
  double RScale; //for TMP profile, Rs/Rs0
  double PotScale; //for TMP profile, Pots/Pots0
  double K; //M=KR, for isothermal profile, where K=Vc^2/G.
  int IsForbidden; //flag telling whether the parameters are forbidden (e.g, negative mass parameter)
//   int TMPid;//for TMP profile
  int virtype;
  HaloType_t type;
} Halo_t;

extern void halo_set_type(HaloType_t t, VirType_t virtype, double Redshift, const double scales[], Halo_t *halo, int TMPid);
extern void halo_set_param(const double pars[], Halo_t *halo);
extern double halo_mass(double r, Halo_t *halo);
extern double halo_pot(double r, Halo_t *halo);

#define HALO_HEADER_INCLUDED
#endif
