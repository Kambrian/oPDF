#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mymath.h"
#include "potential.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>

float OmegaM0=0.3,OmegaL0=0.7,hubble_par0=0.73;

void evolve_cosmology(float z, struct CosmParZ *cosm)
{
//calculate cosmology and virial factors at redshift z, to fill *cosm
float Hratio,scaleF,x;
cosm->z=z;
scaleF=1.0/(1+z);
cosm->Hz=HUBBLE0 * sqrt(OmegaM0 /scaleF/scaleF/scaleF + (1 -OmegaM0 -OmegaL0) / scaleF/scaleF +OmegaL0);
Hratio=cosm->Hz/HUBBLE0;
cosm->OmegaZ=OmegaM0/scaleF/scaleF/scaleF/Hratio/Hratio;
x=cosm->OmegaZ-1;
cosm->virialF[VIR_TH]=18.0*3.1416*3.1416+82.0*x-39.0*x*x;//<Rho_vir>/Rho_cri
cosm->virialF[VIR_C200]=200.;
cosm->virialF[VIR_B200]=200.*cosm->OmegaZ;
}

double comoving_virial_radius(double Mvir, double z, int virtype)
{
//Mvir: in units of 10^10Msun/h
//virial_type being one of VIR_TH, VIR_B200, VIR_C200
//return comoving virial radius in Mpc/h
if(Mvir<0) Mvir=-Mvir; //generalize to negative mass

if(VIR_B200==virtype)	
return cbrt(G*Mvir/OmegaM0/100/HUBBLE0/HUBBLE0);  //R_B200 is independent of z and cosmology

struct CosmParZ cosm;
evolve_cosmology(z,&cosm);
return cbrt(2.0*G*Mvir/cosm.virialF[virtype]/cosm.Hz/cosm.Hz)*(1+z);
}

double NFW_concentration(double M, double z, int virtype)
//concentration parameter for a given mass, redshift and virtype
{
/*	
% M: 10^10Msun/h
% z: redshift
% virtype: virial definition: 0: Bryan-Norman; 1: 200c; 2: 200b.	
* */
double c;
double A,B,C,Mp;
//Duffy et al., 2008, Full Sample. z=0~2.
Mp=2e2; //10^10Msun/h
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

void decode_NFWprof(double z, double M, double c, int virtype, struct NFWParZ *halo)
{//translate halo (M,c) to other profile parameters
 //M: 10^10Msun/h
// (z,virtype) also need to be present upon input.
double scaleF,rhoc;
struct CosmParZ cosm;

halo->z=z;
halo->M=M;
halo->c=c;
halo->virtype=virtype;

evolve_cosmology(halo->z,&cosm);
// scaleF=1.0/(1+halo->z);
rhoc=(3.0*cosm.Hz*cosm.Hz)/(8.0*M_PI*G);
// printf("rhoc=%g\n",rhoc);
halo->Rhos=cosm.virialF[halo->virtype]/3.0*halo->c*halo->c*halo->c/(log(1+halo->c)-halo->c/(1+halo->c))*rhoc;//physical
//halo->Rhos*=scaleF*scaleF*scaleF;//convert to comoving
halo->Rv=cbrt(fabs(halo->M)/(4.0*M_PI/3.0*cosm.virialF[halo->virtype]*rhoc));//physical, generalized to negative mass
//halo->Rv/=scaleF;//convert to comoving
halo->Rs=halo->Rv/halo->c;	//physical
halo->Ms=4.0*M_PI*halo->Rhos*halo->Rs*halo->Rs*halo->Rs;
halo->Pots=-G*halo->Ms/halo->Rs;
// printf("Rhos=%g,Rs=%g,Rv=%g,Pots=%g\n", halo->Rhos, halo->Rs, halo->Rv, halo->Pots);
}

void decode_NFWprof2(double z, double Rhos, double Rs, int virtype, struct NFWParZ *halo)
{//translate halo (Rhos,Rs) to other profile parameters
 //M: 10^10Msun/h
// (z,virtype) also need to be present upon input.
double scaleF,rhoc;
struct CosmParZ cosm;

halo->z=z;
halo->Rhos=Rhos;
halo->Rs=Rs;
halo->virtype=virtype;

evolve_cosmology(halo->z,&cosm);
// scaleF=1.0/(1+halo->z);
rhoc=(3.0*cosm.Hz*cosm.Hz)/(8.0*M_PI*G);
// printf("rhoc=%g\n",rhoc);
halo->Ms=4.0*M_PI*halo->Rhos*halo->Rs*halo->Rs*halo->Rs;
halo->Pots=-G*halo->Ms/halo->Rs;
//to be done:xxxxxxxxxxxxxxxx
halo->M=0.;
halo->c=0.;
halo->Rv=0.;
// printf("Rhos=%g,Rs=%g,Rv=%g,Pots=%g\n", halo->Rhos, halo->Rs, halo->Rv, halo->Pots);
}

void NFWMC_setp(double *pars, Halo_t *halo)
{
#define PropID_Pots 0
#define PropID_Rs   1  
#define PropID_Mass 2
#define PropID_Conc 3
#define PropID_Rhos 4
  memcpy(halo->pars, pars, NumParMax*sizeof(double));
  halo->property[PropID_Mass]=pars[0];
  halo->property[PropID_Conc]=pars[1];
  halo->property[PropID_Rhos]=
#define VIRFactor 200 //200*crit
  double scaleF,rhoc;
  rhoc=(3.0*HUBBLE0*HUBBLE0)/(8.0*M_PI*G);
struct CosmParZ cosm;

evolve_cosmology(halo->z,&cosm);
// scaleF=1.0/(1+halo->z);

// printf("rhoc=%g\n",rhoc);
halo->Rhos=VIRFactor/3.0*halo->property[3]*halo->c*halo->c/(log(1+halo->c)-halo->c/(1+halo->c))*rhoc;//physical
//halo->Rhos*=scaleF*scaleF*scaleF;//convert to comoving
halo->Rv=cbrt(fabs(halo->M)/(4.0*M_PI/3.0*cosm.virialF[halo->virtype]*rhoc));//physical, generalized to negative mass
//halo->Rv/=scaleF;//convert to comoving
halo->Rs=halo->Rv/halo->c;	//physical
halo->Ms=4.0*M_PI*halo->Rhos*halo->Rs*halo->Rs*halo->Rs;
halo->Pots=-G*halo->Ms/halo->Rs;
}
 
double NFW_pot(double r, Halo_t *halo)
{
  
}

void set_halo_type(HaloType_t type, Halo_t *halo)
{
  halo->type=type;
  switch(type)
  {
	case NFWMC:
	  halo->pot=
  }
}

void set_halo_pars(double *pars, Halo_t *halo)
{
  
}