#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mymath.h"
#include "globals.h"
#include "cosmology.h"
#include "halo.h"
#include "template.h"
#include "nfw.h"

void halo_set_type(HaloType_t t, VirType_t virtype, double Redshift, const double scales[], Halo_t *halo, int TMPid)
{
  halo->type=t;
  halo->z=Redshift;
  halo->virtype=virtype;
  memcpy(halo->scales, scales, NUM_PAR_MAX*sizeof(double));
  if(t==HT_TMPMC||t==HT_TMPPotScaleRScale)//use template, allocate it
  {
	int cid=get_current_TMPid();
	if(cid!=TMPid)
	{
	  if(cid>=0)
		fprintf(stderr, "==INFO: Changing halo template from %d to %d\n"
		"Note a maximum of only one template can be loaded at any time.\n"
		"The old template is now replaced\n", cid, TMPid);
	  free_potential_spline();
	  init_potential_spline(TMPid);
	}
  }
}
#define ISPHYSICAL(x) ((x)>0)
void halo_set_param(const double pars[], Halo_t *halo)
{
  memcpy(halo->pars, pars, NUM_PAR_MAX*sizeof(double));
  double phypar[NUM_PAR_MAX];
  int i;
  halo->IsForbidden=0; //You can set the range of allowed parameters by setting halo.IsForbidden=1 when input parameters are not allowed, which will cause the likelihood routines to return -INFINITY for log-like or INFINITY for distances.
  for(i=0;i<NUM_PAR_MAX;i++)
	phypar[i]=pars[i]*halo->scales[i];
  switch(halo->type)
  {
	case HT_NFWMC:
	  halo->M=phypar[0];
	  halo->c=phypar[1];
	  halo->IsForbidden=!(ISPHYSICAL(halo->M)&&ISPHYSICAL(halo->c));
	  if(!halo->IsForbidden) decode_NFWprof(halo);
	  break;
	case HT_NFWPotsRs:
	case HT_CorePotsRs://same as NFWPotsRs
	  halo->Pots=-phypar[0];
	  halo->Rs=phypar[1];
	  halo->Rhos=-halo->Pots/halo->Rs/halo->Rs/4/M_PI/Globals.units.Const.G;
	  halo->Ms=-halo->Pots*halo->Rs/Globals.units.Const.G;
	  halo->IsForbidden=!(ISPHYSICAL(-halo->Pots)&&ISPHYSICAL(halo->Rs));
	  break;
	case HT_NFWRhosRs:
	case HT_CoreRhosRs:
	  halo->Rhos=phypar[0];
	  halo->Rs=phypar[1];
	  halo->Ms=4.0*M_PI*halo->Rhos*halo->Rs*halo->Rs*halo->Rs;
	  halo->Pots=-Globals.units.Const.G*halo->Ms/halo->Rs;
	  halo->IsForbidden=!(ISPHYSICAL(halo->Rhos)&&ISPHYSICAL(halo->Rs));
	  break;
	case HT_TMPMC:
	  halo->M=phypar[0];
	  halo->c=phypar[1];
	  halo->IsForbidden=!(ISPHYSICAL(halo->M)&&ISPHYSICAL(halo->c));
	  if(!halo->IsForbidden) decode_TemplateProf(halo);
	  break;
	case HT_TMPPotScaleRScale:
	  halo->PotScale=phypar[0];
	  halo->RScale=phypar[1];
	  halo->IsForbidden=!(ISPHYSICAL(halo->PotScale)&&ISPHYSICAL(halo->RScale));
	  break;
	default:
	  DEBUGPRINT("Error: Unknown halo parametrization %d\n", halo->type);
	  exit(1);
  }
}

double halo_mass(double r, Halo_t *halo)
{
  if(halo->IsForbidden) return 0.;
  switch(halo->type)
  {
	case HT_TMPMC:
	case HT_TMPPotScaleRScale:
	  return eval_density_spline(r/halo->RScale)*halo->PotScale/halo->RScale/halo->RScale*4.*M_PI/3.*r*r*r; //use spline if inited. rho(r)=A/B^2 rho0(r/B).
	  //bugfix: r/halo->Rs to r/halo->RScale
	case HT_NFWMC:
	case HT_NFWPotsRs:
	case HT_NFWRhosRs:
	  return NFW_mass(r, halo);
	default:
	  DEBUGPRINT("Error: mass profile does not support parametrization %d yet\n", halo->type);
	  exit(1);
  }
}

double halo_pot(double r, Halo_t *halo)
{ 
  if(halo->IsForbidden) return 0.;
  double x=r/halo->Rs;
  switch(halo->type)
  {
	case HT_NFWMC:
	case HT_NFWPotsRs:
	case HT_NFWRhosRs:
	  if(x<EPS) return halo->Pots; //to avoid numerical nan at x=0;
	  return halo->Pots*log(1+x)/x; //the halo boundary is irrelevant in this problem, since only the potential difference affects the orbit
	case HT_CorePotsRs:
	case HT_CoreRhosRs:
	  if(x<EPS) return halo->Pots/2.;
	  return halo->Pots*(log(1.+x)/x-0.5/(1.+x));
	case HT_TMPMC:
	case HT_TMPPotScaleRScale:
	  return eval_potential_spline(r/halo->RScale)*halo->PotScale; 
	default:
	  DEBUGPRINT("Error: halo type %d not supported by halo_pot yet", halo->type);
	  exit(1);
  }
  return 0.;//just to conform with the rules for a final return.
}