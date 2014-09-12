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
 
static int cmpPartE(const void *p1, const void *p2)
{ //in ascending order
  if(((Particle_t *)p1)->E > ((Particle_t *)p2)->E ) 
	return 1;
  
  if(((Particle_t *)p1)->E < ((Particle_t *)p2)->E ) 
	return -1;
  
  return 0;
}
 
static int cmpPartL(const void *p1, const void *p2)
 { //in ascending order
   if(((Particle_t *)p1)->L2 > ((Particle_t *)p2)->L2 ) 
     return 1;
   
   if(((Particle_t *)p1)->L2 < ((Particle_t *)p2)->L2 ) 
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
 
void calibrate_particle_weights(Tracer_t *Sample)
{//calibrate weights, so that sum(weight)=num_particles. also rescale the average mass so that total mass is still correct.
  int i;
  double w;
  
  if(Sample->nP==0) return;
  
  for(i=0,w=0;i<Sample->nP;i++)
	w+=Sample->P[i].w;
  w/=Sample->nP; //average weight
  for(i=0;i<Sample->nP;i++)
	Sample->P[i].w/=w;
  Sample->mP*=w;
}

void load_tracer_particles(char *datafile, Tracer_t * Sample)
{
	size_t nload;
    int i,j, *p;
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
    
	sprintf(B.name,"%sHaloID",grpname);
    nload=load_hdfmatrix(datafile,&B,1,H5T_NATIVE_INT);
	if(nload>0)
	{
	  p=B.x;
	  for(i=0;i<Sample->nP;i++)	Sample->P[i].haloid=p[i];
	  free(B.x);
	}
	
	sprintf(B.name,"%sSubID",grpname);
    nload=load_hdfmatrix(datafile,&B,1,H5T_NATIVE_INT);
	if(nload>0)
	{
	  p=B.x;
	  for(i=0;i<Sample->nP;i++)	Sample->P[i].subid=p[i];
	  free(B.x);
	}
/*	
	sprintf(B.name,"%sStrmID",grpname);
    nload=load_hdfmatrix(datafile,&B,1,H5T_NATIVE_INT);
	if(nload>0)
	{
	  p=B.x;
	  for(i=0;i<Sample->nP;i++)	Sample->P[i].strmid=p[i];
	  free(B.x);
	}
*/    
    sprintf(B.name,"%sflag",grpname);
    nload=load_hdfmatrix(datafile,&B,1,H5T_NATIVE_INT);
    if(nload==0) printf("Assuming flag=1 for every particle\n");
	if(nload>0)
	{
	  p=B.x;
	  for(i=0;i<Sample->nP;i++)	Sample->P[i].flag=p[i];
	  free(B.x);
	}
    else
	  for(i=0;i<Sample->nP;i++)	Sample->P[i].flag=1;
	
	sprintf(A.name,"%sPartMass",grpname);
    nload=load_hdfmatrixF(datafile,&A,1);
	Sample->mP=1.; //initial unit
	if(nload<=1)
	{
	  if(nload==1) 
		Sample->mP=A.x[0];
	  for(i=0;i<Sample->nP;i++)	Sample->P[i].w=1.;
	}
	else if(nload==Sample->nP)
	{
	  for(i=0;i<nload;i++)
		Sample->P[i].w=A.x[i];
	  calibrate_particle_weights(Sample);
	}
	else
	{
	  fprintf(stderr, "Error: incomplete PartMass arr %zd (expecting %d)\n", nload, Sample->nP);
	  exit(1);
	}
	if(nload>0) free(A.x);
	printf("mP=%g\n", Sample->mP);
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

void shuffle_tracer_particles(unsigned long int seed, Tracer_t *Sample)
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

void resample_tracer_particles(unsigned long int seed, Tracer_t *ReSample, Tracer_t *Sample)
{//bootstrap resampling (with replacement)
    const gsl_rng_type * T;
    gsl_rng * r;
	
	copy_tracer_particles(0, -1, ReSample, Sample); //null copy
	ReSample->P=malloc(sizeof(Particle_t)*Sample->nP);
    /* create a generator chosen by the
       environment variable GSL_RNG_TYPE */
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);
    
    gsl_ran_sample(r, ReSample->P, ReSample->nP, Sample->P, Sample->nP, sizeof(Particle_t));
    
    gsl_rng_free (r);
	
	calibrate_particle_weights(ReSample);
}

void copy_tracer_particles(int offset, int sample_size, Tracer_t *Sample, Tracer_t *FullSample)
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
  Sample->FlagUseWeight=FullSample->FlagUseWeight;
  Sample->mP=FullSample->mP;
  calibrate_particle_weights(Sample);
  //its the user's responsibility to manage RadialCount[]
}

void cut_tracer_particles(Tracer_t *Sample, double rmin, double rmax)
{//apply radial cuts, in place
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
  calibrate_particle_weights(Sample);
}

void squeeze_tracer_particles(Tracer_t *Sample)
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
  calibrate_particle_weights(Sample);
}

void free_tracer_particles(Tracer_t *Sample)
{
  if(Sample->nP)
  {
  free(Sample->P);
  Sample->nP=0;
  }
}

void free_tracer(Tracer_t *Sample)
{
  free_tracer_rcounts(Sample);
  free_tracer_views(Sample);
  free_tracer_particles(Sample);
}

void print_tracer_particle(Tracer_t *Sample, int i)
{
  printf("%d, %p\n", Sample->nP, Sample->P);
  printf("%g, %g\n", Sample->P[i].x[0], Sample->P[i].r);
  printf("%g-%g\n", Sample->rmin, Sample->rmax);
}

void sort_part_flag(Particle_t *P, int nP)
{//ascending in flag
  qsort(P, nP, sizeof(Particle_t), cmpPartFlag);
}
void sort_part_L(Particle_t *P, int nP)
{//ascending in L
  qsort(P, nP, sizeof(Particle_t), cmpPartL);
}
void sort_part_E(Particle_t *P, int nP)
{//ascending in E
  qsort(P, nP, sizeof(Particle_t), cmpPartE);
}
void sort_part_R(Particle_t *P, int nP)
{//ascending in R
  qsort(P, nP, sizeof(Particle_t), cmpPartR);
}
void create_tracer_views(Tracer_t *Sample, int nView, char proxy)
{//sort Sample and divide into nView equal-size subsamples, discarding remainders.
  //the data are not copied. so only views are generated.

  if(nView>1)//otherwise no need to sort
  {
	switch(proxy)
	{
	  case 'E':
		qsort(Sample->P, Sample->nP, sizeof(Particle_t), cmpPartE);
		break;
	  case 'L':
		qsort(Sample->P, Sample->nP, sizeof(Particle_t), cmpPartL);
		break;
	  case 'r':
		qsort(Sample->P, Sample->nP, sizeof(Particle_t), cmpPartR);
		break;
	  case 'f':
		qsort(Sample->P, Sample->nP, sizeof(Particle_t), cmpPartFlag);
		break;
	  default:
		fprintf(stderr, "error: unknown proxy=%c\n", proxy);
		exit(1);
	}
  }
  free_tracer_views(Sample); //always try to clean up first.
  Sample->Views=calloc(nView, sizeof(TracerView));//now safe to re-alloc
  int i,nP,offset;
  nP=Sample->nP/nView;
  for(i=0,offset=0;i<nView;i++)
  {
	copy_tracer_particles(0, -1, Sample->Views+i, Sample);
	Sample->Views[i].nP=nP;
	Sample->Views[i].P=Sample->P+offset;
	Sample->Views[i].mP=Sample->mP;
	calibrate_particle_weights(Sample->Views+i);
	if(proxy=='r')
	{
	  Sample->Views[i].rmin=Sample->Views[i].P[0].r;
	  Sample->Views[i].rmax=Sample->Views[i].P[nP-1].r;
	}
	offset+=nP;
  }
  Sample->nView=nView;
  Sample->ViewType=proxy;
}
void free_tracer_views(Tracer_t *Sample)
{
  if(Sample->nView)
  {
	int i;
	for(i=0;i<Sample->nView;i++)//deep cleaning first
	  free_tracer_views(Sample->Views+i);
	free(Sample->Views);
	Sample->nView=0;
	Sample->ViewType='\0';
  }
}
void free_tracer_rcounts(Tracer_t *Sample)
{
  if(Sample->nbin_r)
  {
	free(Sample->RadialCount);
	Sample->nbin_r=0;
  }
}

void count_tracer_radial(Tracer_t *Sample, int nbin, int FlagRLogBin)
{//count inside linear radial bins
  //init rmin-rmax before calling.
  int i,j;
  double dr,logRmin;
  if(Sample->rmax==0) {fprintf(stderr, "Error: cut particles or assign rmin/rmax first\n"); exit(1);}
  if(Sample->nbin_r!=nbin) //need to realloc
	free_tracer_rcounts(Sample);
  Sample->FlagRLogBin=FlagRLogBin;
  Sample->nbin_r=nbin;
  Sample->RadialCount=calloc(Sample->nbin_r, sizeof(double));
  if(Sample->FlagRLogBin)
  {
	logRmin=log(Sample->rmin);
	dr=(log(Sample->rmax)-logRmin)/Sample->nbin_r;
  }
  else
	dr=(Sample->rmax-Sample->rmin)/Sample->nbin_r;
  #pragma omp parallel
  {
	double *bincount;
	bincount=calloc(Sample->nbin_r, sizeof(double));
    #pragma omp for private(j)
    for(i=0;i<Sample->nP;i++)
    { 
	  if(Sample->FlagRLogBin)
		j=floor((log(Sample->P[i].r)-logRmin)/dr);
	  else
		j=floor((Sample->P[i].r-Sample->rmin)/dr);
      if(j<0) j=0;
      if(j>=Sample->nbin_r) j=Sample->nbin_r-1;
	  if(Sample->FlagUseWeight)
		bincount[j]+=Sample->P[i].w;
	  else
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
int SubSampleSize=1000;
void init_tracer(Tracer_t *Sample)
{
  char datafile[1024]=ROOTDIR"/data/mockhalo_wenting.hdf5";
  Sample->rmin=1;
  Sample->rmax=1000;
  HaloM0=183.5017;
  HaloC0=16.1560;
  HaloProfID=0;
  
  if(NULL!=getenv("DynDataFile"))
  {
    printf("Importing parameters from environment..\n");
    sprintf(datafile,"%s/data/%s", ROOTDIR, getenv("DynDataFile"));
    SubSampleSize=strtol(getenv("DynSIZE"),NULL, 10);
    Sample->rmin=strtod(getenv("DynRMIN"),NULL);
    Sample->rmax=strtod(getenv("DynRMAX"),NULL);
    HaloM0=strtod(getenv("DynM0"),NULL);
    HaloC0=strtod(getenv("DynC0"),NULL);
	HaloProfID=strtol(getenv("DynProfID"), NULL, 10);
  }
  else
    printf("Warning: Using default parameters with datafile %s .\n", datafile);
 
  printf("%s; %d; %g,%g;%g,%g\n", datafile, SubSampleSize, Sample->rmin,Sample->rmax,HaloM0,HaloC0);
  decode_NFWprof(HaloZ0,HaloM0,HaloC0,VIR_C200,&Halo);
  HaloRhos0=Halo.Rhos;
  HaloRs0=Halo.Rs;
  
  load_tracer_particles(datafile, Sample);
  cut_tracer_particles(Sample,Sample->rmin,Sample->rmax);
  shuffle_tracer_particles(100,Sample);
  count_tracer_radial(Sample, NumRadialCountBin, 1);
}
void make_sample(int offset, int samplesize, Tracer_t *Sample, Tracer_t *FullSample)
{//make a subsample
  copy_tracer_particles(offset, samplesize, Sample, FullSample);
  count_tracer_radial(Sample, NumRadialCountBin, 1);
}