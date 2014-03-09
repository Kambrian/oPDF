#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mymath.h"
#include "cosmology.h"
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
double distance_modulus(double z)
{//distance modulus, M=m-DM(z)-Kcorr
  return 5.0*log10(lum_dist_flat(z)/hubble_par0/1e-5); //-2.5*log10[(DL/10pc)^2]
}
double lum_dist_flat(double z)
{//luminosity dist for flat universe, 
 //in units of Mpc/h
  return AD_dist_flat(0.,z)*(1.+z)*(1.+z);
}
double AD_dist_flat(double z1,double z2)
{
/*generalized angular diameter distance fitting formula for flat cosmology
% in  units of Mpc/h
% accurate to within 0.3% for OmegaM=(0.2~1) and z=(0.01~inf).
% Ref: M. Adachi and M. Kasai (arXiv:1012.2670; 1111.6396)
% set z1=0 for usual DA
*/

double x1,Fx1,x2,Fx2;

if (1!=OmegaM0+OmegaL0)
{
	printf("error: not flat cosmology for AD_dist_flat\n");
	exit(2);
}

x1=(1.-OmegaM0)/OmegaM0/pow(1.+z1,3);
Fx1=1./sqrt(1.+z1)*(2. + 2.641*x1 + 0.8830*x1*x1 + 0.05313*x1*x1*x1)/(1. + 1.392*x1 + 0.5121*x1*x1 + 0.03944*x1*x1*x1);
x2=(1.-OmegaM0)/OmegaM0/pow(1.+z2,3);
Fx2=1./sqrt(1.+z2)*(2. + 2.641*x2 + 0.8830*x2*x2 + 0.05313*x2*x2*x2)/(1. + 1.392*x2 + 0.5121*x2*x2 + 0.03944*x2*x2*x2);
return 3000./(1.+z2)/sqrt(OmegaM0)*(Fx1-Fx2);
}

double sigma_crit(double zl,double zs)
{
//critical surface density for lensing, physical
double c=3e5; //light speed, km/s
double DL,DS,DLS;
//double chiL,chiS;

if(1==OmegaM0+OmegaL0)
{
DL=AD_dist_flat(0,zl);
DS=AD_dist_flat(0,zs);
DLS=AD_dist_flat(zl,zs);
}
else
{
	printf("error: this sigma_crit function currently only supports flat cosmology\n");
	exit(2);
}
//~ {
//~ chiL=comoving_dist(zl);
//~ chiS=comoving_dist(zs);
//~ DL=comoving_function(chiL,OmgM+OmgL)./(1+zl);
//~ DS=comoving_function(chiS,OmgM+OmgL);  //(1+zs) factor ommitted
//~ DLS=comoving_function(chiS-chiL,OmgM+OmgL); //(1+zs) factor ommitted
//~ }
return c*c/4.0/M_PI/G*DS/DL/DLS;   //10^10Msun/h/(Mpc/h)^2
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
halo->Pots=-4*M_PI*G*halo->Rhos*halo->Rs*halo->Rs;
// printf("Rhos=%g,Rs=%g,Rv=%g,Pots=%g\n", halo->Rhos, halo->Rs, halo->Rv, halo->Pots);
}

void decode_NFWprof2(double z, double Rhos, double Rs, int virtype, struct NFWParZ *halo)
{//translate halo (M,c) to other profile parameters
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
halo->Pots=-4.0*M_PI*G*halo->Rhos*halo->Rs*halo->Rs;
//to be done:xxxxxxxxxxxxxxxx
halo->M=0.;
halo->c=0.;
halo->Rv=0.;
// printf("Rhos=%g,Rs=%g,Rv=%g,Pots=%g\n", halo->Rhos, halo->Rs, halo->Rv, halo->Pots);
}

double NFW_DeltSig(double r,struct NFWParZ *halo)
{
/*	
% DeltaSigma=sigmacrit*shear for NFW halo,physical value, ref: Wright & Brainerd 2000
 //10^10Msun/h/(Mpc/h)^2
% r: Mpc/h, physical
*/

double sig,x,x2;
x=r/halo->Rs;
x2=x*x;

// if(halo->M<0) return -INFINITY;
// if(halo->M==0) return 0.;

#define EPSILON 1e-5
if(fabs(x-1)<EPSILON)
	sig=10.0/3-4*log(2.0);
else if(x<1)
	sig=8*atanh(sqrt((1-x)/(1+x)))/x2/sqrt(1-x2)+4./x2*log(x/2)
        -2./(x2-1)+4*atanh(sqrt((1-x)/(1+x)))/(x2-1)/sqrt(1-x2);
else
	sig=8*atan(sqrt((x-1)/(1+x)))/x2/sqrt(x2-1)+4./x2*log(x/2)
        -2./(x2-1)+4*atan(sqrt((x-1)/(1+x)))/(x2-1)/sqrt(x2-1);   
   
sig*=halo->Rs*halo->Rhos;//physical
#ifdef INCLUDE_2HALO_LENS
double scaleFinv=(1+halo->z);
double D=growth_factor(halo->z);
sig+=halo_bias(halo->M)*lin_DeltSig_comov(r*scaleFinv)*scaleFinv*scaleFinv*D*D; //bias term
#endif
return (halo->M<0?-sig:sig); //generalize to negative mass. this might not be good for MCMC, better return 0 probability in that case, or you need to customize your proposal function.
}
double ISO_DeltSig(double r, struct ISOParZ *halo)
{// to be finished.......................
	return 0.;
}
double halo_DeltSig(double r, struct HaloParZ *halo)
{//physical r input and physical density output
	if(halo->ModelID==100) return 0.; //null model
	
	if(halo->ModelID<10)//NFW
	return NFW_DeltSig(r,&halo->nfw);
	else
	return ISO_DeltSig(r,&halo->iso);
}

double halo_bias(double M)
{/* fitting formula to the linear bias at z=0.2
  input: halo mass M200b, 10^10Msun/h

% accurate to 5% percent
% x=log10(M), y=log10(b)
% biasfit = 
% 
%      Linear model Poly5:
%      biasfit(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
%      Coefficients (with 95% confidence bounds):
%        p1 =   4.378e-06  (3.673e-06, 5.083e-06)
%        p2 =  -3.705e-05  (-6.717e-05, -6.94e-06)
%        p3 =  -0.0003215  (-0.000785, 0.000142)
%        p4 =    0.003687  (0.0005968, 0.006777)
%        p5 =    -0.02313  (-0.0315, -0.01476)
%        p6 =    -0.07584  (-0.08276, -0.06892)
*/

if(M<0) return 0.;
double x=log10(M)+10.; //change to Msun/h  
if(x>17.) x=17.;
if(x<0.) x=0.;
double b=POLY5(4.378e-6, -3.705e-5, -0.0003215, 0.003687, -0.02313, -0.07584, x);
return pow(10.,b);
}
double lin_DeltSig_comov(double r)
{/* input: comoving r, Mpc/h
  * output: linearDsig, comoving, 10^10Msun/h/(Mpc/h)^2, for z=0 and b=1.
  % well fitted by polynomial, within 5 percent (0.02dex).
% y=log10(s), x=log10(r):
% corrfit = 
% 
%      Linear model Poly5:
%      corrfit(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
%      Coefficients (with 95% confidence bounds):
%        p1 =    -0.00325  (-0.003798, -0.002701)
%        p2 =    -0.02437  (-0.02591, -0.02283)
%        p3 =    -0.05605  (-0.05912, -0.05298)
%        p4 =      -0.156  (-0.1622, -0.1498)
%        p5 =      0.3879  (0.3823, 0.3935)
%        p6 =       1.518  (1.513, 1.522) */
double x=log10(r);
double s=POLY5(-0.00325,-0.02437,-0.05605,-0.156,0.3879,1.518,x);
return pow(10.,s);
}
double growth_factor(double z)
{/* return growth factor normalized to z=0; applicable for z<1.
  % x=z, y=D.
% well fitted to within 0.1 percent (extremely accurate)
% growthfit = 
% 
%      Linear model Poly2:
%      growthfit(x) = p1*x^2 + p2*x + p3
%      Coefficients (with 95% confidence bounds):
%        p1 =      0.1318  (0.1313, 0.1324)
%        p2 =     -0.5195  (-0.52, -0.5189)
%        p3 =           1  (1, 1)*/
return POLY2(0.1318,-0.5195,1.,z);
}