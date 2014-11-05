#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// #include "mymath.h"
#include "globals.h"
#include "cosmology.h"
#include "halo.h"
#include "template.h"
#include "nfw.h"

void halo_set_type(HaloType_t t, Halo_t *halo)
{
  halo->type=t;
}
void halo_set_param(double *pars, Halo_t *halo)
{
  halo->z=Globals.cosmology.Redshift;
  halo->virtype=Globals.virtype;
//   memcpy(Halo.pars, pars, NUM_PAR_MAX*sizeof(double));
  switch(halo->type)
  {
	case HT_NFWMC:
	  halo->M=pars[0];
	  halo->c=pars[1];
	  decode_NFWprof(halo);
	  break;
	case HT_NFWPotsRs:
	case HT_CorePotsRs://same as NFWPotsRs
	  halo->Pots=-pars[0];
	  halo->Rs=pars[1];
	  halo->Rhos=-halo->Pots/halo->Rs/halo->Rs/4/M_PI/Globals.units.Const.G;
	  halo->Ms=-halo->Pots*halo->Rs/Globals.units.Const.G;
	  break;
	case HT_NFWRhosRs:
	case HT_CoreRhosRs:
	  halo->Rhos=pars[0];
	  halo->Rs=pars[1];
	  halo->Ms=4.0*M_PI*halo->Rhos*halo->Rs*halo->Rs*halo->Rs;
	  halo->Pots=-Globals.units.Const.G*halo->Ms/halo->Rs;
	  break;
	case HT_TMPMC:
	  halo->M=pars[0];
	  halo->c=pars[1];
	  decode_TemplateProf(halo);
	case HT_TMPPotScaleRScale:
	  halo->PotScale=pars[0];
	  halo->RScale=pars[1];
	  break;
	default:
	  DEBUGPRINT("Error: Unknown halo parametrization %d\n", halo->type);
	  exit(1);
  }
}

double halo_mass(double r, Halo_t *halo)
{
  switch(halo->type)
  {
	case HT_TMPMC:
	case HT_TMPPotScaleRScale:
	  return eval_density_spline(r/halo->Rs)*halo->PotScale/halo->RScale/halo->RScale*4.*M_PI/3.*r*r*r; //use spline if inited.
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