#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mymath.h"
#include "cosmology.h"

double NFW_concentration(double M, double z, int virtype)
//concentration parameter for a given mass, redshift and virtype
{
/*	
% M: virial mass
% z: redshift
% virtype: virial definition: 0: Bryan-Norman; 1: 200c; 2: 200b.	
* */
double c;
double A,B,C,Mp;
//Duffy et al., 2008, Full Sample. z=0~2.
Mp=2e12/Units.MassInMsunh; //10^10Msun/h
switch(virtype)
{
	case VIR_TH:
	A=7.85;
	B=-0.081;
	C=-0.71;
	break;
	case VIR_C200:
	A=5.71;
	B=-0.084;
	C=-0.47;
	break;
	case VIR_B200:
	A=10.14;
	B=-0.081;
	C=-1.01;
	break;
	default:
	fprintf(stderr,"error: unknown virial type %d\n", virtype);
	exit(1);
}

if(M<0) M=-M; //generalize to negative mass

c=A*pow(M/Mp,B)*pow(1+z,C);

return c;
}

void decode_NFWprof(Halo_t *halo)
{//translate halo (M,c) to other profile parameters
// (z,virtype) also need to be present upon input.
double scaleF,rhoc;
struct CosmParZ cosm;

evolve_cosmology(halo->z,&cosm);
// scaleF=1.0/(1+halo->z);
rhoc=(3.0*cosm.Hz*cosm.Hz)/(8.0*M_PI*Units.Const.G);
// printf("rhoc=%g\n",rhoc);
halo->Rhos=cosm.virialF[halo->virtype]/3.0*halo->c*halo->c*halo->c/(log(1+halo->c)-halo->c/(1+halo->c))*rhoc;//physical
//halo->Rhos*=scaleF*scaleF*scaleF;//convert to comoving
halo->Rv=cbrt(fabs(halo->M)/(4.0*M_PI/3.0*cosm.virialF[halo->virtype]*rhoc));//physical, generalized to negative mass
//halo->Rv/=scaleF;//convert to comoving
halo->Rs=halo->Rv/halo->c;	//physical
halo->Ms=4.0*M_PI*halo->Rhos*halo->Rs*halo->Rs*halo->Rs;
halo->Pots=-Units.Const.G*halo->Ms/halo->Rs;
// printf("Rhos=%g,Rs=%g,Rv=%g,Pots=%g\n", halo->Rhos, halo->Rs, halo->Rv, halo->Pots);
}

double NFW_mass(double r, Halo_t *halo)
{ 
  double x=r/halo->Rs;
  return halo->Ms*(log(1+x)-x/(1+x)); 
}