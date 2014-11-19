#include <math.h>
#include <stdio.h>
#include <gsl/gsl_spline.h>

#include "mymath.h"
#include "globals.h"
#include "cosmology.h"
#include "halo.h"
#include "template.h"

struct SplineData
{
	int FlagUseSpline; //spline ready to be used. override the default potential calculation with spline interpolation
	gsl_interp_accel *acc;
	gsl_spline *spline ;
	gsl_interp_accel *acc_dens;
	gsl_spline *spline_dens;
	int TMPid;
	double Rs;
};
static struct SplineData PotSpline;
#pragma omp threadprivate(PotSpline)  //make sure each thread has its own cache
int get_current_TMPid()
{
  if(PotSpline.FlagUseSpline)
	return PotSpline.TMPid;
  
  return -1;
}
void init_potential_spline(int TMPid)
{
#include "TemplateData.h"
	if(TMPid<0)
	{
	  DEBUGPRINT("Error: TMPid=%d, no profile data\n", TMPid);
	  exit(1);
	}
	#pragma omp parallel 
	{
		PotSpline.acc= gsl_interp_accel_alloc ();
		PotSpline.spline= gsl_spline_alloc (gsl_interp_cspline, LEN_PROF);
		gsl_spline_init(PotSpline.spline, PotentialTemplate[TMPid][0], PotentialTemplate[TMPid][1], LEN_PROF);
		
		PotSpline.acc_dens= gsl_interp_accel_alloc ();
		PotSpline.spline_dens= gsl_spline_alloc (gsl_interp_cspline, LEN_PROF);
		gsl_spline_init(PotSpline.spline_dens, PotentialTemplate[TMPid][0], PotentialTemplate[TMPid][2], LEN_PROF);
		
		PotSpline.FlagUseSpline=1;
		PotSpline.TMPid=TMPid;
		PotSpline.Rs=TemplateScale[TMPid];
	}
// 	printf("Spline: [%g,%g], %zd, x~[%g, %g], y~[%g,%g]\n", PotSpline.spline->interp->xmin, PotSpline.spline->interp->xmax,
// 			PotSpline.spline->size, PotSpline.spline->x[0], PotSpline.spline->x[PotSpline.spline->size-1], 
// 									PotSpline.spline->y[0], PotSpline.spline->y[PotSpline.spline->size-1]);
}
void free_potential_spline()
{
	#pragma omp parallel
	{
	 if(PotSpline.FlagUseSpline)
	 {
	  gsl_spline_free (PotSpline.spline);
	  gsl_interp_accel_free(PotSpline.acc);
	  
	  gsl_spline_free (PotSpline.spline_dens);
	  gsl_interp_accel_free(PotSpline.acc_dens);
	  
	  PotSpline.TMPid=-1;
	  PotSpline.FlagUseSpline=0;
	 }
	}
}
double eval_potential_spline(double r)
{
  if(r>PotSpline.spline->interp->xmax)//outside the interpolation range: zero density, potential from the inner sphere
  {
	int imax=PotSpline.spline->size-1;
	return -PotSpline.spline->y[imax]*PotSpline.spline->x[imax]/r;
  }
  return -gsl_spline_eval(PotSpline.spline, r, PotSpline.acc); //the input are |pot|, restore the sign.
}

double eval_density_spline(double r)
{//cumulative density
  if(r>PotSpline.spline_dens->interp->xmax)//outside the interpolation range: zero density, cumdensity from the inner sphere	
  {
	int imax=PotSpline.spline_dens->size-1;
	double x=PotSpline.spline_dens->x[imax]/r;
	return PotSpline.spline_dens->y[imax]*x*x*x;
  }
  return gsl_spline_eval(PotSpline.spline_dens, r, PotSpline.acc_dens);
}

void decode_TemplateProf(Halo_t *halo)
{//translate halo (M,c) to other profile parameters (PotScale, RScale). Input are physical values.
double scaleF,rhoc;
struct CosmParZ cosm;

evolve_cosmology(halo->z,&cosm);
rhoc=(3.0*cosm.Hz*cosm.Hz)/(8.0*M_PI*Globals.units.Const.G);
halo->Rv=cbrt(fabs(halo->M)/(4.0*M_PI/3.0*cosm.virialF[halo->virtype]*rhoc));
halo->Rs=halo->Rv/halo->c; 
halo->RScale=halo->Rs/PotSpline.Rs;
halo->PotScale=cosm.virialF[halo->virtype]*rhoc*halo->RScale*halo->RScale/eval_density_spline(halo->Rv/halo->RScale);
}