//to test the numerical convergence of TS calculation, as a function of TOLs
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
  int estimator=10,subsample_id=0,i; 
  if(argc!=3)
  {
    fprintf(stderr, "Usage: %s [estimator] [sampleid]\nNow exit.\n",argv[0]);
    exit(1);
  }
  estimator=atoi(argv[1]);
  subsample_id=atoi(argv[2]);
  
/*  double x,en=sqrt(1000.);
  for(x=0.022;x<0.023;x+=0.0001)
  {
    printf("%f,%g:%f\n",x,(en+0.12+0.11/en)*x,KSprob(1000,x));
  }
  return 0;
*/   
  init();
  select_particles(subsample_id);
//   printf("%g,%g\n", MODEL_TOL_BIN, MODEL_TOL_REL);
//   freeze_and_like(pars, estimator);
  for(MODEL_TOL_BIN=1e-1;MODEL_TOL_BIN>1e-5;MODEL_TOL_BIN/=10)
  {
	  printf("%g:\t", MODEL_TOL_BIN);
	  for(MODEL_TOL_REL=1e-3;MODEL_TOL_REL>1e-6;MODEL_TOL_REL/=10)
	  {printf("%g, %g; ", MODEL_TOL_REL, freeze_and_like(pars, estimator));fflush(stdout);}
	  printf("\n");
  }
/*
  time_t t1,t2;
  t1=time(NULL);
  MODEL_TOL_BIN=1e-2;MODEL_TOL_REL=1e-2;
  for(pars[0]=0.95;pars[0]<1.05;pars[0]+=0.01)
  {printf("%g,", freeze_and_like(pars,estimator));fflush(stdout);}
  t2=time(NULL);
  printf("\n %ld sec\n", t2-t1);fflush(stdout);
  t1=time(NULL);
  MODEL_TOL_BIN=1e-2;MODEL_TOL_REL=1e-3;//this is pretty good
  for(pars[0]=0.95;pars[0]<1.05;pars[0]+=0.01)
  {printf("%g,", freeze_and_like(pars,estimator));fflush(stdout);}
  t2=time(NULL);
  printf("\n %ld sec\n", t2-t1);fflush(stdout);
  t1=time(NULL);
  MODEL_TOL_BIN=1e-2;MODEL_TOL_REL=1e-4;
  for(pars[0]=0.95;pars[0]<1.05;pars[0]+=0.01)
  {printf("%g,", freeze_and_like(pars,estimator));fflush(stdout);}
  t2=time(NULL);
  printf("\n %ld sec\n", t2-t1);fflush(stdout);
  t1=time(NULL);
  MODEL_TOL_BIN=1e-3;MODEL_TOL_REL=1e-3;
  for(pars[0]=0.95;pars[0]<1.05;pars[0]+=0.01)
  {printf("%g,", freeze_and_like(pars,estimator));fflush(stdout);}
  t2=time(NULL);
  printf("\n %ld sec\n", t2-t1);fflush(stdout);
  t1=time(NULL);
  MODEL_TOL_BIN=1e-4;MODEL_TOL_REL=1e-7;
  for(pars[0]=0.95;pars[0]<1.05;pars[0]+=0.01)
  { printf("%g,", freeze_and_like(pars,estimator));fflush(stdout);}
  t2=time(NULL);
  printf("\n %ld sec\n", t2-t1);fflush(stdout);
  t1=time(NULL);*/
  free_integration_space();
  
  return 0;
  
}