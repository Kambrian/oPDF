#ifndef HALO_HEADER_INCLUDED

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
} HaloType_t;

typedef struct Halo
{//all quantities are physical
//   double pars[NUM_PAR_MAX]; //raw parameters
  double z;
  double M;
  double c;
  double Rv;
  double Pots;//-4*pi*G*rhos*rs^2, the potential at r=0
  double Rs;
  double Rhos;//4*pi*rs^3*rhos
  double Ms;
  double RScale; //for TMP profile, Rs/Rs0
  double PotScale; //for TMP profile, Pots/Pots0
  int TMPid;//for TMP profile
  int virtype;
  HaloType_t type;
} Halo_t;

extern void halo_set_type(HaloType_t t, Halo_t *halo);
extern void halo_set_param(double *pars, Halo_t *halo);
extern double halo_mass(double r, Halo_t *halo);
extern double halo_pot(double r, Halo_t *halo);

#define HALO_HEADER_INCLUDED
#endif
