
extern NFWParZ Halo;
#ifndef POTENTIAL_HEADER_INCLUDED
/*constants in alternative GADGET units
 * Mass: 10^10Msun *
 * Length: kpc *
 * vel:   km/s *  */
#define G 43007.1
#define HUBBLE0 0.073  //h=0.73, H0=h*0.1 (km/s/kpc)

/* virial definitions */
#define VIR_TH 0
#define VIR_C200 1
#define VIR_B200 2

#define NumPropMax 10

typedef enum
{
  NFWMC=0,
  NFWPotsRs,
  NFWRhosRs,
  GNFWCusp,
} HaloType_t;

struct Halo;
typedef struct Halo Halo_t;
typedef void HaloSetFunc_t(double *, Halo_t *);
typedef double PotFunc_t(double, Halo_t *);
typedef double DensityFunc_t(double, Halo_t *);
typedef double MassProfFunc_t(double, Halo_t *);
struct Halo
{
  double property[NumPropMax]; //Halo properties converted from pars
  HaloType_t type;
  HaloSetFunc_t * setp;
  PotFunc_t * pot;
  DensityFunc_t * density;
  MassProfFunc_t * mass;
};


#define POTENTIAL_HEADER_INCLUDED
#endif
