#include <math.h>

#include "cosmology.h"

void evolve_cosmology(float z, struct CosmParZ *cosm)
{
//calculate cosmology and virial factors at redshift z, to fill *cosm
float Hratio,scaleF,x;
cosm->z=z;
scaleF=1.0/(1+z);
cosm->Hz=Globals.units.Const.H0 * sqrt(Globals.cosmology.OmegaM0 /scaleF/scaleF/scaleF + (1 -Globals.cosmology.OmegaM0 -Globals.cosmology.OmegaL0) / scaleF/scaleF +Globals.cosmology.OmegaL0);
Hratio=cosm->Hz/Globals.units.Const.H0;
cosm->OmegaZ=Globals.cosmology.OmegaM0/scaleF/scaleF/scaleF/Hratio/Hratio;
x=cosm->OmegaZ-1;
cosm->virialF[VIR_TH]=18.0*3.1416*3.1416+82.0*x-39.0*x*x;//<Rho_vir>/Rho_cri
cosm->virialF[VIR_C200]=200.;
cosm->virialF[VIR_B200]=200.*cosm->OmegaZ;
}

double comoving_virial_radius(double Mvir, double z, VirType_t virtype)
{
//virial_type being one of VIR_TH, VIR_B200, VIR_C200
//return comoving virial radius
if(Mvir<0) Mvir=-Mvir; //generalize to negative mass
if(VIR_B200==virtype)	
return cbrt(Globals.units.Const.G*Mvir/Globals.cosmology.OmegaM0/100./Globals.units.Const.H0/Globals.units.Const.H0);  //R_B200 is independent of z and cosmology

struct CosmParZ cosm;
evolve_cosmology(z,&cosm);
return cbrt(2.0*Globals.units.Const.G*Mvir/cosm.virialF[virtype]/cosm.Hz/cosm.Hz)*(1.+z);
}
