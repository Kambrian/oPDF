//create mock satellites from DM mimicking sat orbits
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

int SAMPLE_SIZE_FACTOR;
#define DATADIR "/cosma/home/durham/jvbq85/data/SatDyn/data"

extern void mock_stars(int HaloId, int seed, Tracer_t *NewStar);
extern void save_mockstars(int HaloId, Tracer_t *NewStar, hid_t file_id);
int main(int argc, char **argv)
{
  float hubble=0.73;
  set_units(1e10*hubble, hubble, 1);

  SAMPLE_SIZE_FACTOR=1;
  if(argc>2) SAMPLE_SIZE_FACTOR=atoi(argv[2]);
  
  hid_t file_id;
  char datafile[1024];
  sprintf(datafile,"%s/IsolatedMWHalosSublist.FromDM.Factor%d.hdf5", DATADIR, SAMPLE_SIZE_FACTOR);
  if(try_readfile(datafile))
      file_id=H5Fopen(datafile, H5F_ACC_RDWR, H5P_DEFAULT);
  else
      file_id = H5Fcreate (datafile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file

  int HaloId;
 // for(HaloId=0;HaloId<944;HaloId++)
  {
    HaloId=atoi(argv[1]);
      printf("Halo %d==============================\n", HaloId);
    Tracer_t NewStar={};
    mock_stars(HaloId, 100, &NewStar);
    save_mockstars(HaloId, &NewStar, file_id);
    free_tracer(&NewStar);
  }
  
  H5Fclose(file_id);
  
  return 0;
}

void read_halo_param(int HaloId, double pars[2])
{
    char datafile[1024];
    sprintf(datafile, "%s/for_Zhaozhou/haloparams_isolated.txt", DATADIR);
    FILE *fp;
    myfopen(fp, datafile, "r");
    char ignore[1024];
    fgets(ignore, 1024, fp);//skip header
    int i;
    for(i=0;i<HaloId;i++)
        fgets(ignore, 1024, fp);
    double mass_err;
    fscanf(fp, "%lf %lf %lf", pars, &mass_err, pars+1);
    printf("%d, %f, %f\n", HaloId, pars[0], pars[1]);
    fclose(fp);
}

int get_MaxStarInEView(Tracer_t *Star, int istart, TracerView *View)
{
    int i;
    for(i=istart;i<Star->nP;i++)
        if(Star->P[i].E>View->proxybin[1]) break;
    return i;
}

int get_MaxStarInLView(Tracer_t *Star, int istart, int iend, TracerView *View)
{
    int i;
    for(i=istart;i<iend;i++)
        if(Star->P[i].L2>View->proxybin[1]) break;
    return i;
}

void mock_stars(int HaloId, int seed, Tracer_t *NewStar)
{
  int AddHubbleFlow=0;
  Tracer_t DM={}, Star={};
  NewStar->nP=0; NewStar->nView=0;
  char datafile[1024];
  //double rmin=1.,rmax=400.;
   
  sprintf(datafile,"%s/IsolatedMWHalos.3Rvir.hdf5",DATADIR);
  char grpname[256];
  sprintf(grpname, "Halo%03d/", HaloId);
  load_tracer_particles(datafile, grpname, &DM, AddHubbleFlow);
  //cut_tracer_particles(&DM,rmin,rmax);
//   shuffle_tracer_particles(100, &DM);
//   DM.nP=10000;

  sprintf(datafile,"%s/IsolatedMWHalosSublist.hdf5", DATADIR);
  load_tracer_particles(datafile, grpname, &Star, AddHubbleFlow);
  //cut_tracer_particles(&Star,rmin,rmax);
//   shuffle_tracer_particles(100, &Star);
//   Star.nP=1000;
  int i;
  for(i=0;i<DM.nP;i++)
	DM.P[i].flag=(DM.P[i].subid<=0);
  squeeze_tracer_particles(&DM);
  for(i=0;i<Star.nP;i++)
	Star.P[i].flag=(Star.P[i].subid<=0);
  squeeze_tracer_particles(&Star);
  printf("DM: %d, Star: %d\n", DM.nP, Star.nP);
  
  NewStar->nP=Star.nP*SAMPLE_SIZE_FACTOR;
  NewStar->P=malloc(sizeof(Particle_t)*NewStar->nP);
  
  double scales[NUM_PAR_MAX]={100.,1.}, pars[NUM_PAR_MAX]={1.,1.};
  read_halo_param(HaloId, pars);
  //halo_set_type(HT_TMPPotScaleRScale, VIR_C200, 0., scales, &(DM.halo), TMPid);
  halo_set_type(HT_NFWMC, VIR_C200, 0., scales, &(DM.halo), 0);
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
  int pid0=0, eid;
  NewStar->nP=0;
  for(eid=0; eid<nbin[0]; eid++)
  {
//       printf("\neid=%d\t",eid);
      int pid;
      pid=get_MaxStarInEView(&Star, pid0, &DM.Views[eid]);
      if(eid==nbin[0]-1) pid=Star.nP;//last bin right open [ )
      if(pid>pid0)//open
      {
          sort_part_L(Star.P+pid0, pid-pid0);
          int cid0=pid0, lid;
          for(lid=0;lid<nbin[1];lid++)
          {
              int cid;
              TracerView *View=&(DM.Views[eid].Views[lid]);
              cid=get_MaxStarInLView(&Star, cid0, pid, View);
              if(lid==nbin[1]-1) cid=pid;
              if(cid>cid0)//sample
              {
//                   printf("%d:%d,%d\t", lid, cid0, cid);
                  int nP=cid-cid0;
                  nP*=SAMPLE_SIZE_FACTOR;
                  if(nP>View->nP)
                  {
                    printf("Warning: nstar=%d, nDM=%d\n", nP, View->nP); 
                    nP=View->nP;
                  }
                  gsl_ran_choose(r, NewStar->P+NewStar->nP, nP, View->P, View->nP, sizeof(Particle_t));
                  NewStar->nP+=nP;
                  cid0=cid;
              }
          }
          pid0=pid;
      }
  }
  gsl_rng_free (r);
  printf("Old: %d; New: %d\n", Star.nP, NewStar->nP);
  NewStar->P=realloc(NewStar->P, sizeof(Particle_t)*NewStar->nP);
  
  free_tracer(&DM);
  free_tracer(&Star);
}
void save_mockstars(int HaloId, Tracer_t *NewStar, hid_t file_id)
{
  hid_t grp_id;
  herr_t status;
  hsize_t dims[2];
  /*char datafile[1024];
  sprintf(datafile,"%s/IsolatedMWHalosGalaxylist.FromDM.hdf5", DATADIR);
  if(try_readfile(datafile))
      file_id=H5Fopen(datafile, H5F_ACC_RDWR, H5P_DEFAULT);
  else
      file_id = H5Fcreate (datafile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file
*/
  char grpname[1024];
  sprintf(grpname, "Halo%03d", HaloId);
  //H5Gunlink(file_id, grpname);
  grp_id=H5Gcreate(file_id, grpname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
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
  status = H5LTmake_dataset(grp_id,"x",2,dims,H5T_NATIVE_FLOAT,pos);
  status = H5LTmake_dataset(grp_id,"v",2,dims,H5T_NATIVE_FLOAT,vel);
  dims[1]=1;
  status = H5LTmake_dataset(grp_id,"SubID",2,dims,H5T_NATIVE_INT,subid);
  free(pos);
  free(vel);
  free(subid);
  
  H5Gclose(grp_id);
 // H5Fclose(file_id);
}
