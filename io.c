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
//FIXME: the calibration of weights in views destroy the relative weighting of particles accross views! so do not do weighted likelihood once you created views.
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
	else if(proxy=='E')
	{
	  Sample->Views[i].proxybin[0]=Sample->Views[i].P[0].E;
	  if(i==nView-1)
		Sample->Views[i].proxybin[1]=Sample->Views[i].P[nP-1].E;
	  else
		Sample->Views[i].proxybin[1]=Sample->Views[i].P[nP].E; //first point in next bin; so Ebin[0]<=E<Ebin[1].
	}
	else if(proxy=='L')
	{
	  Sample->Views[i].proxybin[0]=Sample->Views[i].P[0].L2;
	  if(i==nView-1)
		Sample->Views[i].proxybin[1]=Sample->Views[i].P[nP-1].L2;
	  else
		Sample->Views[i].proxybin[1]=Sample->Views[i].P[nP].L2; //first point in next bin; so Ebin[0]<=E<Ebin[1].
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
void mock_stars(char halo, int seed, Tracer_t *NewStar)
{
  Tracer_t DM={}, Star={};
  NewStar->nP=0; NewStar->nView=0;
  char datafile[1024];
  double rmin=10.,rmax=300.;
  switch(halo)
  {
	case 'A':
	  HaloM0=184.2;
	  HaloC0=16.10;
	  HaloProfID=2;
	  break;
	case 'B':
	  HaloM0=81.94;
	  HaloC0=8.16;
	  HaloProfID=3;
	  break;
	case 'C':
	  HaloM0=177.4;
	  HaloC0=12.34;
	  HaloProfID=4;
	  break;
	case 'D':
	  HaloM0=177.4;
	  HaloC0=8.73;
	  HaloProfID=5;
	  break;
	case 'E':
	  HaloM0=118.5;
	  HaloC0=8.67;
	  HaloProfID=6;
	  break;
	default:
	  printf("Error: wrong halo %c\n",halo);
	  exit(1);
  }
  decode_NFWprof(HaloZ0,HaloM0,HaloC0,VIR_C200,&Halo);
  HaloRhos0=Halo.Rhos;
  HaloRs0=Halo.Rs;
 
  sprintf(datafile,"%s/data/FullData/%c2DM.hdf5", ROOTDIR, halo);
  load_tracer_particles(datafile, &DM);
  cut_tracer_particles(&DM,rmin,rmax);
//   shuffle_tracer_particles(100, &DM);
//   DM.nP=10000;

  sprintf(datafile,"%s/data/%c2star.hdf5", ROOTDIR, halo);
  load_tracer_particles(datafile, &Star);
  cut_tracer_particles(&Star,rmin,rmax);
//   shuffle_tracer_particles(100, &Star);
//   Star.nP=1000;
  int i;
  for(i=0;i<DM.nP;i++)
	DM.P[i].flag=(DM.P[i].subid<=0);
  squeeze_tracer_particles(&DM);
  for(i=0;i<Star.nP;i++)
	Star.P[i].flag=(Star.P[i].subid<=0);
  squeeze_tracer_particles(&Star);
  
  NewStar->nP=Star.nP;
  NewStar->P=malloc(sizeof(Particle_t)*NewStar->nP);
  
  init_potential_spline();
  double pars[2]={1.,1.};
  int nbin[2]={100,100};
  char proxies[32]="EL";
  create_nested_views(pars, nbin, proxies, &DM);
  freeze_energy(pars, &Star);
  free_potential_spline();

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);
	
  qsort(Star.P, Star.nP, sizeof(Particle_t), cmpPartE);
  int pid, pid0, eid;
  NewStar->nP=0;
  for(eid=0,pid0=0,pid=0;pid<Star.nP&&eid<nbin[0];pid++)
  {
	while(Star.P[pid].E>=DM.Views[eid].proxybin[1]||pid==Star.nP-1)//passed bin boundary, process and walk eid
	{
// 	  printf("\n%d:",eid);
	  if(pid>pid0)//non-empty bin, open it
	  {
		qsort(Star.P+pid0, pid-pid0, sizeof(Particle_t), cmpPartL);
		int cid,cid0, lid;
		for(lid=0,cid0=pid0,cid=pid0;(cid<pid&&lid<nbin[1])||cid==Star.nP-1;cid++)
		{
		  Tracer_t *View=&(DM.Views[eid].Views[lid]);
		  while(Star.P[cid].L2>=View->proxybin[1]||cid==Star.nP-1)//passed a new bin, sample and walk lid
		  {
// 			printf("%d,",lid);
			int nP=cid-cid0;
			if(nP>0)//non-empty, sample it
			{
			  if(nP>View->nP)
			  {
				printf("Warning: nstar=%d, nDM=%d\n", nP, View->nP); 
				nP=View->nP;
			  }
			  gsl_ran_choose(r, NewStar->P+NewStar->nP, nP, View->P, View->nP, sizeof(Particle_t));
			  NewStar->nP+=nP;
			  cid0=cid;
			}
			lid++;
			if(lid==nbin[1]) break;
			View=&(DM.Views[eid].Views[lid]);
		  }
		}
		pid0=pid;//rebase lower boundary
	  }
	  eid++;
	  if(eid==nbin[0]) break;
	}
  }
  gsl_rng_free (r);
  printf("Old: %d; New: %d\n", Star.nP, NewStar->nP);
  NewStar->P=realloc(NewStar->P, sizeof(Particle_t)*NewStar->nP);
  
  free_tracer(&DM);
  free_tracer(&Star);
}
void save_mockstars(char halo, Tracer_t *NewStar)
{
  hid_t file_id;
  herr_t status;
  hsize_t dims[2];
  char datafile[1024];
  sprintf(datafile,"%s/data/%c2starCleanFromDM.hdf5", ROOTDIR, halo);
  file_id = H5Fcreate (datafile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file

  dims[0]=NewStar->nP;
  dims[1]=3;
  float (* pos)[3], (* vel)[3];
  int * subid;
  int i,j;
  pos=malloc(sizeof(float)*3*NewStar->nP);
  vel=malloc(sizeof(float)*3*NewStar->nP);
  subid=malloc(sizeof(int)*NewStar->nP);
  for(i=0;i<NewStar->nP;i++)
  {
	for(j=0;j<3;j++)
	{
	  pos[i][j]=NewStar->P[i].x[j];
	  vel[i][j]=NewStar->P[i].v[j];
	}
	subid[i]=NewStar->P[i].subid;
  }
  status = H5LTmake_dataset(file_id,"/x",2,dims,H5T_NATIVE_FLOAT,pos);
  status = H5LTmake_dataset(file_id,"/v",2,dims,H5T_NATIVE_FLOAT,vel);
  dims[1]=1;
  status = H5LTmake_dataset(file_id,"/SubID",2,dims,H5T_NATIVE_INT,subid);
  free(pos);
  free(vel);
  free(subid);
  
  H5Fclose(file_id);
}