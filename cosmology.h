#ifndef COSMOLOGY_HEADER_INCLUDED

#include "globals.h"

struct CosmParZ
{
	double z;
	double OmegaZ;
	double Hz;
	double virialF[3];
};

/*cosmology params*/
extern void evolve_cosmology(float z, struct CosmParZ *cosm);
extern double comoving_virial_radius(double Mvir, double z, VirType_t virial_type);
#define COSMOLOGY_HEADER_INCLUDED
#endif
