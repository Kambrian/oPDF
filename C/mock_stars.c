#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/times.h>

#include "mymath.h"
#include "cosmology.h"
#include "globals.h"
#include "halo.h"
#include "tracer.h"
#include "models.h"
#include "hdf_util.h"

#define DATADIR "/gpfs/data/wenting/data"

extern void mock_stars(char *HaloName, int TMPid, int seed, Tracer_t *NewStar);
extern void save_mockstars(char *HaloName, Tracer_t *NewStar);
int main(int argc, char **argv)
{
  Tracer_t NewStar={};
  if(argc!=3)
  {
	printf("Usage: %s [HaloName] [TemplateID]\nNow exit.\n", argv[0]);
	exit(1);
  }
  char HaloName[1024];
  strcpy(HaloName, argv[1]);
  int TMPid=atoi(argv[2]);
  
  mock_stars(HaloName, TMPid, 100, &NewStar);
  save_mockstars(HaloName, &NewStar);
  free_tracer(&NewStar);
  return 0;
}

void mock_stars(char *HaloName, int TMPid, int seed, Tracer_t *NewStar)
{
  int AddHubbleFlow=0;
  Tracer_t DM={}, Star={};
  NewStar->nP=0; NewStar->nView=0;
  char datafile[1024];
  double rmin=10.,rmax=300.;
   
  sprintf(datafile,"%s/%sDM.hdf5",DATADIR, HaloName);
  load_tracer_particles(datafile, &DM, AddHubbleFlow);
  cut_tracer_particles(&DM,rmin,rmax);
//   shuffle_tracer_particles(100, &DM);
//   DM.nP=10000;

  sprintf(datafile,"%s/%sStar.hdf5", DATADIR, HaloName);
  load_tracer_particles(datafile, &Star, AddHubbleFlow);
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
  
  double scales[NUM_PAR_MAX]={1.,1.}, pars[NUM_PAR_MAX]={1.,1.};
  halo_set_type(HT_TMPPotScaleRScale, VIR_C200, 0., scales, &(DM.halo), TMPid);
  halo_set_param(pars, &(DM.halo));
  Star.halo=DM.halo;
  tracer_set_energy(&DM);
  tracer_set_energy(&Star);
  int nbin[2]={100,100};
  char proxies[32]="EL";
  create_nested_views(nbin, proxies, &DM);
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);
	
  sort_part_E(Star.P, Star.nP);
  int pid, pid0, eid;
  NewStar->nP=0;
  for(eid=0,pid0=0,pid=0;pid<Star.nP&&eid<nbin[0];pid++)
  {
	while(Star.P[pid].E>=DM.Views[eid].proxybin[1]||pid==Star.nP-1)//passed bin boundary, process and walk eid
	{
// 	  printf("\n%d:",eid);
	  if(pid>pid0)//non-empty bin, open it
	  {
		sort_part_L(Star.P+pid0, pid-pid0);
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
void save_mockstars(char *HaloName, Tracer_t *NewStar)
{
  hid_t file_id;
  herr_t status;
  hsize_t dims[2];
  char datafile[1024];
  sprintf(datafile,"%s/%sStarCleanFromDM.hdf5", DATADIR, HaloName);
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