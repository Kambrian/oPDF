#ifndef COSMOLOGY_HEADER_INCLUDED
/*constants in alternative GADGET units
 * Mass: 10^10Msun *
 * Length: kpc *
 * vel:   km/s *  */
#define G 43007.1
#define HUBBLE0 0.073  //h=0.73, H0=h*0.1 (km/s/kpc)

#define NUM_PAR_MAX 10

/* virial definitions */
#define VIR_TH 0
#define VIR_C200 1
#define VIR_B200 2

struct CosmParZ
{
	double z;
	double OmegaZ;
	double Hz;
	double virialF[3];
};

struct NFWParZ
{
	double pars[NUM_PAR_MAX]; //raw parameters
	double z;
	double M; //Mass, 10^10Msun/h
	double c;
	double Rv; //physical virial radius, Mpc/h
	double Rs; //physical scale radius, Mpc/h
	double Rhos; //physical scale density
	double Pots; //-4*pi*G*rhos*rs^2, the potential at r=0
	double Ms; //4*pi*rs^3*rhos
	int virtype;
};
struct ISOParZ //isothermal parameter
{
	double z;
	double M;
	double vdisp;
	double Rv;
	double Rhos;
	int virtype;
};
struct HaloParZ
{
	int ModelID;//0~10: NFW; >10: ISO.
	struct NFWParZ nfw;
	struct ISOParZ iso;
};
/*cosmology params*/
extern float OmegaM0,OmegaL0,hubble_par0;

extern void evolve_cosmology(float z, struct CosmParZ *cosm);
extern double comoving_virial_radius(double Mvir, double z, int virial_type);
extern double z_selection_frac(double z);
extern double distance_modulus(double z);
extern double lum_dist_flat(double z);
extern double AD_dist_flat(double z1,double z2);
extern double sigma_crit(double zl,double zs);
extern double NFW_concentration(double M, double z, int virtype);
extern void decode_NFWprof(double z, double M, double c, int virtype, struct NFWParZ *halo);
extern void decode_NFWprof2(double z, double Rhos, double Rs, int virtype, struct NFWParZ *halo);
extern double NFW_DeltSig(double r,struct NFWParZ *halo);
extern double ISO_DeltSig(double r, struct ISOParZ *halo);
extern double halo_DeltSig(double r, struct HaloParZ *halo);
extern double halo_bias(double M);
extern double lin_DeltSig_comov(double r);
extern double growth_factor(double z);
#define COSMOLOGY_HEADER_INCLUDED
#endif
