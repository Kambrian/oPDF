//tailor wenting's data
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mymath.h"
#include "hdf_util.h"

#define h0 0.73
#define LUNIT ((1/h0)*1000) //kpc, for B4, A2, B2,....

struct PList
{
  int np;
  int *PIndex;
};
extern void load_particles(char halo);
extern void free_particles();
extern void collect_particles(struct PList * p);
extern void dump_particles_hdf(char *outfile, struct PList *p);
float *CoM, * VCoM;
int main(int argc, char** argv)
{
	char outfile[1024];

	struct PList p;
	float centers[5][3]={{57.06108093,  52.61642456,  48.70379639},
	{ 49.64647675,  50.13929367,  49.77927399},
	{ 50.29379654,  50.32624054,  50.51562119},
	{ 60.87985229,  54.10048676,  52.26412964},
	{ 46.0791626,   76.08740997,  66.36212158}};
	float halovels[5][3]={{ 294.80819702,  100.80763245,  -55.95608902},
	{-306.80792236,  -52.49833298,  -24.17752457},
	{ 340.02838135,  403.45211792,  289.01052856},
	{ 517.9788208,   162.09895325,  128.54304504},
	{ 209.11463928, -170.64677429,  274.96737671}};
	char halo_names[6]="ABCDE";
	
	if(argc!=2)
	{printf("usage: %s [halo_id=0~4]\n",argv[0]);fflush(stdout);exit(1);}
	
	int halo_id=atoi(argv[1]);
	char halo=halo_names[halo_id];
	printf("%d,%c\n",halo_id,halo);
	CoM=centers[halo_id];
	VCoM=halovels[halo_id];
	
	load_particles(halo);

	sprintf(outfile,"/gpfs/data/jvbq85/DynDistr/data/%c2star.hdf5",halo);
	collect_particles(&p);
	dump_particles_hdf(outfile,&p);
	free(p.PIndex);

	free_particles();
	return 0;
}
struct 
{
float (* Pos)[3], (*Vel)[3];
int nP;
} Particles;
void load_particles(char halo)
{
  FloatMat A[2];
  GenericMat B;
  char file[1024];
  sprintf(file, "/gpfs/data/wenting/starpop_output/APCBower/Aq_%c2_WMAP1_Cold/tags_han.hdf5", halo);
  
  printf("loading %s...\n", file);fflush(stdout);
  strcpy(A[0].name,"/PartType0/Coordinates");
  strcpy(A[1].name,"/PartType0/Velocities");
  load_hdfmatrixF(file,A,2);
  if(A[0].size[1]!=3)
  {
	printf("Error, unexpected matrix size %d,%d\n",(int)A[0].size[0],(int)A[0].size[1]);
	exit(1);
  }
  Particles.nP=A[0].size[0];
  Particles.Pos=(float (*)[3])A[0].x;
  Particles.Vel=(float (*)[3])A[1].x;
  printf("%g,%g,%g\n", Particles.Pos[10][0], Particles.Pos[10][1], Particles.Pos[10][2]);
}
void free_particles()
{
free(Particles.Pos);
free(Particles.Vel);
}
	
void dump_particles_hdf(char *outfile, struct PList * p)
{
    hid_t    file_id,group_id;
    herr_t      status;
    size_t     i, j,nread=0;
    hsize_t dims[2];

    //if(try_readfile(outfile))
    file_id = H5Fcreate (outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file
//     group_id = HDFcreate_group(file_id, "/PartType0");
//     status = H5Gclose(group_id);

    float (* pos)[3],(* vel)[3];
    pos=malloc(sizeof(float)*3*p->np);
    vel=malloc(sizeof(float)*3*p->np);
    for(i=0;i<p->np;i++)
    {
      for(j=0;j<3;j++)
      {
	pos[i][j]=(Particles.Pos[p->PIndex[i]][j]-CoM[j])*LUNIT;
	vel[i][j]=Particles.Vel[p->PIndex[i]][j]-VCoM[j];
      }
    }
    dims[0]=p->np;
    dims[1]=3;
    status = H5LTmake_dataset(file_id,"/x",2,dims,H5T_NATIVE_FLOAT,pos);
    status = H5LTmake_dataset(file_id,"/v",2,dims,H5T_NATIVE_FLOAT,vel);
	H5LTset_attribute_float(file_id, "/x", "x0", CoM, 3);
	H5LTset_attribute_float(file_id, "/v", "v0", VCoM, 3);
	
	dims[0]=p->np;
    dims[1]=1;
    status = H5LTmake_dataset(file_id,"/Pid",2,dims,H5T_NATIVE_INT,p->PIndex);
		
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[1]: (%g,%g,%g)\n",pos[1][0],pos[1][1],pos[1][2]);
    free(pos);
    free(vel);
}

typedef float HBTReal;
typedef int HBTInt;
typedef HBTReal HBTxyz[3];
#define BOXSIZE 100.0
typedef HBTReal * AccessPosFunc(HBTInt pid,void *PosData);
typedef struct 
{
	HBTInt ndiv;
	HBTInt np;
	HBTInt UseFullBox; //whether to use a local enclosing-box or full simulation box (latter for periodic boudary)
	HBTInt *hoc;
	HBTInt *list;
	HBTReal range[3][2];
	HBTReal step[3];
	void *PosData; //to be passed to GetPos() to access particle positions
	AccessPosFunc *GetPos;
} LINKLIST;
/*-------linkedlists-------------*/
HBTReal * GetArrPos(HBTInt pid,void *data)
{//pass Pdat.Pos to data
	return ((HBTxyz *) data)[pid];
}
HBTInt linklist_round_gridid(HBTInt i,HBTInt ndiv)
{//to correct for rounding error near boundary 
	if(i<0) return 0;
	if(i>=ndiv) return ndiv-1;
	return i;
}
HBTInt linklist_shift_gridid(HBTInt i,HBTInt ndiv)
{//to correct for periodic conditions; 
//only applicable when def PERIODIC_BDR and ll.UseFullBox=1
	i=i%ndiv;
	if(i<0) i+=ndiv;
	return i;
}
HBTInt linklist_fix_gridid(HBTInt i, LINKLIST *ll)
{
	#ifdef PERIODIC_BDR
	if(ll->UseFullBox)
	return linklist_shift_gridid(i,ll->ndiv);	
	#endif
	return linklist_round_gridid(i,ll->ndiv);
}
HBTInt linklist_get_hoc(LINKLIST *ll, HBTInt i,HBTInt j,HBTInt k)
{
	return ll->hoc[i+j*ll->ndiv+k*ll->ndiv*ll->ndiv];
}
HBTInt linklist_get_hoc_safe(LINKLIST *ll, HBTInt i,HBTInt j,HBTInt k)
{//force fixing of gridids
	#define FIXGRID(i) linklist_fix_gridid(i,ll)
	return ll->hoc[FIXGRID(i)+FIXGRID(j)*ll->ndiv+FIXGRID(k)*ll->ndiv*ll->ndiv];
}

void make_linklist(LINKLIST *ll, HBTInt np,HBTInt ndiv, void *PosData, 
								AccessPosFunc *GetPos, HBTInt UseFullBox)
{
	HBTInt i,j,grid[3];
	HBTInt ind,ndiv2;
	HBTReal x;
	//~ float range[3][2],step[3];
	printf("creating linked list..\n");
	
	ndiv2=ndiv*ndiv;
	ll->ndiv=ndiv;
	ll->np=np;
	ll->UseFullBox=UseFullBox;
	ll->hoc=malloc(sizeof(HBTInt)*ndiv*ndiv*ndiv);
	ll->list=malloc(sizeof(HBTInt)*np);
	ll->PosData=PosData;
	ll->GetPos=GetPos;
	/*determining enclosing cube*/
	if(UseFullBox)
	{
		for(i=0;i<3;i++)
		{
			ll->range[i][0]=0.;
			ll->range[i][1]=BOXSIZE;
		}
		for(j=0;j<3;j++)
			ll->step[j]=BOXSIZE/ndiv;	
	}
	else
	{
		for(i=0;i<3;i++)
			for(j=0;j<2;j++)
				ll->range[i][j]=GetPos(0,PosData)[i];
		for(i=1;i<np;i++)
			for(j=0;j<3;j++)
			{
				x=GetPos(i,PosData)[j];
				if(x<ll->range[j][0])
					ll->range[j][0]=x;
				else if(x>ll->range[j][1])
					ll->range[j][1]=x;
			}
		for(j=0;j<3;j++)
			ll->step[j]=(ll->range[j][1]-ll->range[j][0])/ll->ndiv;
	}
	/*initialize hoc*/
	HBTInt *phoc=ll->hoc;
	for(i=0;i<ndiv*ndiv*ndiv;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((GetPos(i,PosData)[j]-ll->range[j][0])/ll->step[j]);
			grid[j]=linklist_fix_gridid(grid[j],ll);
		}
		ind=grid[0]+grid[1]*ndiv+grid[2]*ndiv2;
		ll->list[i]=ll->hoc[ind];
		ll->hoc[ind]=i; /*use hoc[ind] as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
}

void free_linklist(LINKLIST *ll)
{
	free(ll->hoc);
	free(ll->list);
	ll->ndiv=0;
	ll->np=0;
	ll->PosData=NULL;
	ll->GetPos=NULL;
}

HBTInt * linklist_search_sphere(LINKLIST *ll, HBTReal radius, HBTReal searchcenter[3], HBTInt *N_guess_and_found)
{
	HBTReal dr;
	HBTInt i,j,k,pid,*PIDfound,nfound,nmax;
	int subbox_grid[3][2];
	
	nmax=*N_guess_and_found;
	PIDfound=malloc(sizeof(HBTInt)*nmax);
	nfound=0;
		
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((searchcenter[i]-radius-ll->range[i][0])/ll->step[i]);
	subbox_grid[i][1]=floor((searchcenter[i]+radius-ll->range[i][0])/ll->step[i]);
#ifndef PERIODIC_BDR //do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
	subbox_grid[i][0]=linklist_fix_gridid(subbox_grid[i][0],ll);
	subbox_grid[i][1]=linklist_fix_gridid(subbox_grid[i][1],ll);
#endif	
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc_safe(ll,i,j,k); //in case the grid-id is out of box, in the periodic case
// 				pid=linklist_get_hoc(ll,i,j,k);
				while(pid>=0)
				{
					dr=distance(ll->GetPos(pid,ll->PosData),searchcenter);
					if(dr<radius)
					{
						if(nfound==nmax)
						{
							nmax*=2;
							PIDfound=realloc(PIDfound,sizeof(HBTInt)*nmax);
						}
						PIDfound[nfound]=pid;					
						nfound++;
					}
					pid=ll->list[pid];
				}
			}
	*N_guess_and_found=nfound;		
	return PIDfound;		
}

void collect_particles(struct PList * p)
{
  LINKLIST ll;
  make_linklist(&ll, Particles.nP, 50, Particles.Pos, GetArrPos, 0);
  p->np=Particles.nP;
  p->PIndex=linklist_search_sphere(&ll, 500./LUNIT, CoM, &p->np);
  printf("Found %d particles\n", p->np);
  free_linklist(&ll);
}