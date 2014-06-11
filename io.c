#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "mymath.h"
#include "hdf_util.h"
#include "io.h"
#include "models.h"

static int cmpPartFlag(const void *p1, const void *p2)
 { //in ascending order
   if(((Particle_t *)p1)->flag > ((Particle_t *)p2)->flag ) 
     return 1;
   
   if(((Particle_t *)p1)->flag < ((Particle_t *)p2)->flag ) 
     return -1;
   
   return 0;
 }
 
 static int cmpPartR(const void *p1, const void *p2)
 { //in ascending order
   if(((Particle_t *)p1)->r > ((Particle_t *)p2)->r ) 
     return 1;
   
   if(((Particle_t *)p1)->r < ((Particle_t *)p2)->r ) 
     return -1;
   
   return 0;
 }
void load_tracer(char *datafile, Tracer_t * Sample)
{
    int i,j;
    FloatMat A;
    GenericMat B;
    char grpname[32];
    int dataid=-1;
    if(dataid<0)
      sprintf(grpname,"/");
    else
      sprintf(grpname,"/sample%d/",dataid);
    
    printf("loading %s...\n",datafile);
    sprintf(A.name,"%sx",grpname);
    load_hdfmatrixF(datafile,&A,1);
    if(A.size[1]!=3)
    {
      printf("Error, unexpected matrix size %d,%d\n",(int)A.size[0],(int)A.size[1]);
      exit(1);
    }
	Sample->nP=A.size[0];
    Sample->P=malloc(sizeof(Particle_t)*Sample->nP);
    for(i=0;i<Sample->nP;i++)
    {
      for(j=0;j<3;j++)
      Sample->P[i].x[j]=A.x[i*3+j];
    }
    free(A.x);
    
    sprintf(A.name,"%sv",grpname);
    load_hdfmatrixF(datafile,&A,1);
    if(A.size[1]!=3||A.size[0]!=Sample->nP)
    {
      printf("Error, unexpected matrix size %zd,%zd (expecting %d,3)\n",A.size[0],A.size[1],Sample->nP);
      exit(1);
    }
    for(i=0;i<Sample->nP;i++)
    {
      for(j=0;j<3;j++)
      Sample->P[i].v[j]=A.x[i*3+j];
    }
    free(A.x);
    
    
    sprintf(B.name,"%sflag",grpname);
    size_t nload;
    nload=load_hdfmatrix(datafile,&B,1,H5T_NATIVE_INT);
    if(nload==0) printf("Assuming flag=1 for every particle\n");
    int *p=B.x;
    for(i=0;i<Sample->nP;i++)
      if(nload>0)
	Sample->P[i].flag=p[i];
      else
	Sample->P[i].flag=1;
      
      /*    
    double x0[3]={0.},v0[3]={0.};
    for(i=0;i<Sample->nP;i++)
    {
      for(j=0;j<3;j++)
      {
	x0[j]+=Sample->P[i].x[j];
	v0[j]+=Sample->P[i].v[j];
      }
    }
    for(j=0;j<3;j++)
    {
      x0[j]/=Sample->nP;
      v0[j]/=Sample->nP;
    }
    for(i=0;i<Sample->nP;i++)//shift to center
    {
      for(j=0;j<3;j++)
      {
	Sample->P[i].x[j]-=x0[j];
	Sample->P[i].v[j]-=v0[j];
      }
    }
*/    
    #define VecNorm(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    #define VecProd(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
    for(i=0;i<Sample->nP;i++)
    {
      Sample->P[i].r=sqrt(VecNorm(Sample->P[i].x));
      Sample->P[i].K=VecNorm(Sample->P[i].v)/2.;
      Sample->P[i].L2=cross_product_norm2(Sample->P[i].x,Sample->P[i].v);
      Sample->P[i].vr=VecProd(Sample->P[i].x,Sample->P[i].v)/Sample->P[i].r; //radial vel
    }
    //sort
  // qsort(Sample->P, Sample->nP, sizeof(Particle_t), cmpPartR); //sort according to r
    printf("%d particles loaded\n",Sample->nP);
}

void shuffle_tracer(unsigned long int seed, Tracer_t *Sample)
{
    const gsl_rng_type * T;
    gsl_rng * r;

    /* create a generator chosen by the
       environment variable GSL_RNG_TYPE */
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed+3);
    
    gsl_ran_shuffle(r, Sample->P, Sample->nP, sizeof(Particle_t));
    
    gsl_rng_free (r);
}

void resample_tracer(unsigned long int seed, Tracer_t *ReSample, Tracer_t *Sample)
{//bootstrap resampling (with replacement)
    const gsl_rng_type * T;
    gsl_rng * r;
	
	copy_tracer(0, -1, ReSample, Sample); //null copy
	ReSample->P=malloc(sizeof(Particle_t)*Sample->nP);
    /* create a generator chosen by the
       environment variable GSL_RNG_TYPE */
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);
    
    gsl_ran_sample(r, ReSample->P, ReSample->nP, Sample->P, Sample->nP, sizeof(Particle_t));
    
    gsl_rng_free (r);
}

void copy_tracer(int offset, int sample_size, Tracer_t *Sample, Tracer_t *FullSample)
{
  if(sample_size<0) //null copy
  {
	Sample->nP=0;
	Sample->P=NULL;
  }
  else
  {
	if(sample_size==0) //full copy
	  Sample->nP=FullSample->nP;
	else //sub copy
	  Sample->nP=sample_size;

	if(offset+Sample->nP>FullSample->nP) 
	{
	  printf("error: subsample overflow;%d~%d,nP=%d\n", offset, offset+Sample->nP, FullSample->nP); 
	  exit(1);
	}

	Sample->P=malloc(sizeof(Particle_t)*Sample->nP);
	memcpy(Sample->P,FullSample->P+offset,sizeof(Particle_t)*Sample->nP);
  }
  Sample->rmin=FullSample->rmin;
  Sample->rmax=FullSample->rmax;
  Sample->nbin_r=0;  //means RadialCount un-initialized
}

void cut_tracer(Tracer_t *Sample, double rmin, double rmax)
{//apply radial cuts
  int i,j;
  Sample->rmin=rmin;
  Sample->rmax=rmax;
  for(i=0,j=0;i<Sample->nP;i++)
  {
    if(Sample->P[i].r>rmin&&Sample->P[i].r<rmax)
    {
      if(i>j) Sample->P[j]=Sample->P[i];
      j++;
    }
  }
  Sample->nP=j;
  Sample->P=realloc(Sample->P,sizeof(Particle_t)*Sample->nP);
}

void squeeze_tracer(Tracer_t *Sample)
{//remove flag=0 particles from P, in place.
  int i,j;
  
  for(i=0,j=0;i<Sample->nP;i++)
  {
    if(Sample->P[i].flag)
    {
      if(i>j) Sample->P[j]=Sample->P[i];
      j++;
    }
  }
  Sample->nP=j;
  Sample->P=realloc(Sample->P,sizeof(Particle_t)*Sample->nP);
}

void free_tracer(Tracer_t *Sample)
{
  if(Sample->nbin_r)
  {
	free(Sample->RadialCount);
	Sample->nbin_r=0;
  }
  if(Sample->nP)
  {
  free(Sample->P);
  Sample->nP=0;
  }
}

// void free_data()
// {
//   free_tracer(&DefaultSample);
//   free_tracer(&FullSample);
//   free_integration_space();
// }
void print_tracer(Tracer_t *Sample)
{
  printf("%d, %p\n", Sample->nP, Sample->P);
  printf("%g, %g\n", Sample->P[0].x[0], Sample->P[0].r);
  printf("%g, %g\n", Sample->P[1].x[0], Sample->P[1].r);
  printf("%g-%g\n", Sample->rmin, Sample->rmax);
}

void sort_tracer_flag(Tracer_t *Sample)
{//sort according to flag in ascending order
  qsort(Sample->P, Sample->nP, sizeof(Particle_t), cmpPartFlag); 
}

void count_tracer_radial(Tracer_t *Sample, int nbin)
{//count inside linear radial bins
  int i,j;
  double dr;
  Sample->nbin_r=nbin;
  dr=(Sample->rmax-Sample->rmin)/Sample->nbin_r;
  Sample->RadialCount=calloc(Sample->nbin_r, sizeof(int));
  #pragma omp parallel
  {
	int *bincount;
	bincount=calloc(Sample->nbin_r, sizeof(int));
    #pragma omp for private(j)
    for(i=0;i<Sample->nP;i++)
    {  
      j=floor((Sample->P[i].r-Sample->rmin)/dr);
      if(j<0) j=0;
      if(j>=Sample->nbin_r) j=Sample->nbin_r-1;
      bincount[j]++;
    }
    #pragma omp critical  //note this is not barriered by default!!
    {
    for(j=0;j<Sample->nbin_r;j++)
      Sample->RadialCount[j]+=bincount[j];
    }
    free(bincount);
  }
}
//======higher level funcs====================
int NumRadialCountBin=30;
int SUBSAMPLE_SIZE=1000;
void init_tracer(Tracer_t *Sample)
{
  char datafile[1024]=ROOTDIR"/data/mockhalo_wenting.hdf5";
  Sample->rmin=1;
  Sample->rmax=1000;
  HaloM0=183.5017;
  HaloC0=16.1560;
  
  if(NULL!=getenv("DynDataFile"))
  {
    printf("Importing parameters from environment..\n");
    sprintf(datafile,"%s/data/%s", ROOTDIR, getenv("DynDataFile"));
    SUBSAMPLE_SIZE=strtol(getenv("DynSIZE"),NULL, 10);
    Sample->rmin=strtod(getenv("DynRMIN"),NULL);
    Sample->rmax=strtod(getenv("DynRMAX"),NULL);
    HaloM0=strtod(getenv("DynM0"),NULL);
    HaloC0=strtod(getenv("DynC0"),NULL);
  }
  else
    printf("Warning: Using default parameters with datafile %s .\n", datafile);
 
  printf("%s; %d; %g,%g;%g,%g\n", datafile, SUBSAMPLE_SIZE, Sample->rmin,Sample->rmax,HaloM0,HaloC0);
  decode_NFWprof(HaloZ0,HaloM0,HaloC0,VIR_C200,&Halo);
  HaloRhos0=Halo.Rhos;
  HaloRs0=Halo.Rs;
  
  load_tracer(datafile, Sample);
  cut_tracer(Sample,Sample->rmin,Sample->rmax);
  shuffle_tracer(100,Sample);
  count_tracer_radial(Sample, NumRadialCountBin);
}
void make_sample(int sample_id, Tracer_t *Sample, Tracer_t *FullSample)
{//make a subsample
  copy_tracer(sample_id*SUBSAMPLE_SIZE,SUBSAMPLE_SIZE, Sample, FullSample);
  count_tracer_radial(Sample, NumRadialCountBin);
}