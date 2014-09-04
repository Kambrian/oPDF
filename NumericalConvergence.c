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
// #include "wenting.h"

void test_scan(double tol_bin, double tol_rel, int estimator, Tracer_t *Sample)
{
  time_t t1,t2;
  double pars[NUM_PAR_MAX]={1.,1.};
  t1=time(NULL);
  MODEL_TOL_BIN=tol_bin;MODEL_TOL_REL=tol_rel;
  printf("TOL=%g,%g:\n", tol_bin, tol_rel);
  for(pars[0]=0.95;pars[0]<1.05;pars[0]+=0.01)
  { printf("%g,", freeze_and_like(pars,estimator,Sample));}
  t2=time(NULL);
  printf("\n %ld sec\n", t2-t1);fflush(stdout);
}
int main(int argc, char **argv)
{
  FILE *fp;
  char buf[1024],outdir[1024];
  double pars[NUM_PAR_MAX]={1.,1.,1,1,1,1};
  int estimator=10,samplesize=10000,i; 
  if(argc!=3)
  {
    fprintf(stderr, "Usage: %s [estimator] [samplesize]\nNow exit.\n",argv[0]);
    exit(1);
  }
  estimator=atoi(argv[1]);
  samplesize=atoi(argv[2]);
   
  Tracer_t FullSample={}, Sample={}; //initialize to zeros!!!
  init_tracer(&FullSample);
  make_sample(0, samplesize, &Sample, &FullSample);
  free_tracer(&FullSample);
  
//   dataprob_model(1,1,1, pars);
//   wenting_like(pars, &Sample);
  
  alloc_integration_space();
  
//   printf("%g,%g\n", MODEL_TOL_BIN, MODEL_TOL_REL);
//   freeze_and_like(pars, estimator);
//   MODEL_TOL_BIN_ABS=1e-10;
//   choose_integral_routines(1,0,1);
//   for(MODEL_TOL_BIN=1e-2;MODEL_TOL_BIN>1e-6;MODEL_TOL_BIN/=10)
//   {
// 	  printf("%g:\t", MODEL_TOL_BIN);
// 	  for(MODEL_TOL_REL=1e-3;MODEL_TOL_REL>1e-6;MODEL_TOL_REL/=10)
// 	  {printf("%g, %g; ", MODEL_TOL_REL, freeze_and_like(pars, estimator,&Sample));fflush(stdout);}
// 	  printf("\n");
//   }
//   puts("************************\n");
//     choose_integral_routines(2,0,1);
//     for(MODEL_TOL_BIN=1e-2;MODEL_TOL_BIN>1e-6;MODEL_TOL_BIN/=10)
//   {
// 	  printf("%g:\t", MODEL_TOL_BIN);
// 	  for(MODEL_TOL_REL=1e-3;MODEL_TOL_REL>1e-6;MODEL_TOL_REL/=10)
// 	  {printf("%g, %g; ", MODEL_TOL_REL, freeze_and_like(pars, estimator,&Sample));fflush(stdout);}
// 	  printf("\n");
//   }
//   puts("************************\n");
//   choose_integral_routines(2,0,0);
//     for(MODEL_TOL_BIN=1e-2;MODEL_TOL_BIN>1e-6;MODEL_TOL_BIN/=10)
//   {
// 	  printf("%g:\t", MODEL_TOL_BIN);
// 	  for(MODEL_TOL_REL=1e-3;MODEL_TOL_REL>1e-6;MODEL_TOL_REL/=10)
// 	  {printf("%g, %g; ", MODEL_TOL_REL, freeze_and_like(pars, estimator,&Sample));fflush(stdout);}
// 	  printf("\n");
//   }
//   puts("************************\n");
//     choose_integral_routines(2,1,0);
//   for(MODEL_TOL_BIN=1e-2;MODEL_TOL_BIN>1e-6;MODEL_TOL_BIN/=10)
//   {
// 	  printf("%g:\t", MODEL_TOL_BIN);
// 	  for(MODEL_TOL_REL=1e-3;MODEL_TOL_REL>1e-6;MODEL_TOL_REL/=10)
// 	  {printf("%g, %g; ", MODEL_TOL_REL, freeze_and_like(pars, estimator,&Sample));fflush(stdout);}
// 	  printf("\n");
//   }
/* */
#define scan(a,b) test_scan(a, b, estimator, &Sample)
// choose_integral_routines(1,1,1);scan(1e-4,1e-4);scan(1e-5,1e-7);//slowest, and wrong?
// choose_integral_routines(1,0,1);scan(1e-5,1e-6);scan(1e-5,1e-7);//also wrong
// choose_integral_routines(1,1,0);scan(1e-4,1e-4);scan(1e-5,1e-7);
// choose_integral_routines(1,0,0);scan(1e-4,1e-4);scan(1e-5,1e-7);
// puts("--------------------\n");
// choose_integral_routines(2,1,1);scan(1e-2,1e-3);scan(1e-4,1e-6);
// choose_integral_routines(2,0,1);scan(1e-5,1e-3);//scan(1e-4,1e-6);//dscard.
// choose_integral_routines(2,1,0);
// scan(1e-3,1e-3);
// scan(1e-5,1e-3);//scan(1e-4,1e-6);
// scan(1e-6,1e-6);
// puts("-------------\n");
// scan(1e-2,1e-3);
// scan(1e-3,1e-3);
// scan(1e-4,1e-3);
// scan(1e-5,1e-3);
scan(1e-6,1e-2);
scan(1e-6,1e-3); //speed independent on TOL_BIN
puts("--\n");
scan(1e-6,1e-4);
scan(1e-6,1e-5);
scan(1e-6,1e-6);
//   scan(1e-2,1e-2);
//   scan(1e-2,1e-3);
//   scan(1e-2,1e-4);
//   scan(1e-3,1e-4);
//   scan(1e-4,1e-4);
//   scan(1e-6,1e-8);

  free_integration_space();
  free_tracer(&Sample);
// 	free_potential_spline();
  return 0;
  
}
