#ifndef COSMOLOGY_HEADER_INCLUDED
/*constants in alternative GADGET units
 * Mass: 10^10Msun *
 * Length: kpc *
 * vel:   km/s *  */

/* virial definitions */
typedef enum
{ VIR_TH=0,
  VIR_C200,
  VIR_B200
} VirType_t;

struct CosmParZ
{
	double z;
	double OmegaZ;
	double Hz;
	double virialF[3];
};

/*cosmology params*/
extern float OmegaM0,OmegaL0,RedShift;
extern int VirType;
extern Units_t Units;
extern void set_cosmology(float omegaM0, float omegaL0);
extern void set_units(float MassInMsunh, float LengthInKpch, float VelInKms);
extern void evolve_cosmology(float z, struct CosmParZ *cosm);
extern double comoving_virial_radius(double Mvir, double z, int virial_type);
#define COSMOLOGY_HEADER_INCLUDED
#endif
