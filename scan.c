#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <sys/times.h>

#include "mymath.h"
#include "cosmology.h"
#include "io.h"
#include "models.h"

int main(int argc, char **argv)
{
  FILE *fp;
  char buf[1024],outdir[1024];
  double pars[NUM_PAR_MAX]={1.,1.};
  int subsample_id=0,i;  
  if(argc>=2)
    subsample_id=atoi(argv[1]);
  
/*  double x,en=sqrt(1000.);
  for(x=0.022;x<0.023;x+=0.0001)
  {
    printf("%f,%g:%f\n",x,(en+0.12+0.11/en)*x,KSprob(1000,x));
  }
  return 0;
*/  
  
  double x,y,xlim[3]={-3,3,0.04};//{lims, step}
  sprintf(outdir,"/gpfs/data/jvbq85/DynDistr/data/scanHuge%d",subsample_id);
//   double x,y,xlim[3]={-0.2,0.2,0.002};//{lims, step}
//   sprintf(outdir,"/gpfs/data/jvbq85/DynDistr/data/scanZoom%d",subsample_id);
  mkdir(outdir,0755);
  sprintf(buf,"%s/model%d.T%d",outdir,ESTIMATOR,PHASE_PERIOD);
  myfopen(fp,buf,"w");
  printf("%s (mid=%d): saving to %s...\n",argv[0],ESTIMATOR,buf);fflush(stdout);
  
  init(-1);
  select_particles(subsample_id);
  
  
  for(x=xlim[0];x<=xlim[1];x+=xlim[2])
  {
    for(y=xlim[0];y<xlim[1];y+=xlim[2])
    {
      pars[0]=pow(10.,x);
      pars[1]=pow(10.,y);
      fprintf(fp,"%g\t%g\t%g\n",pars[0],pars[1],freeze_and_like(pars));
    }
    fflush(fp);
  }
  fclose(fp);

  free_integration_space();
  
  return 0;
  
}