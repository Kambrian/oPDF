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
  Tracer_t NewStar={};
  char halo=argv[1][0];
  mock_stars(halo, 100, &NewStar);
  save_mockstars(halo, &NewStar);
  free_tracer(&NewStar);
//   free_tracer(&FullSample);
  return 0;
  
}