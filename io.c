#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mymath.h"
#include "hdf_util.h"
#include "io.h"
  
static Particle *Pall;
static int nPall;
Particle *P;
int nP;

 static int cmpPartR(const void *p1, const void *p2)
 { //in ascending order
   if(((Particle *)p1)->r > ((Particle *)p2)->r ) 
     return 1;
   
   if(((Particle *)p1)->r < ((Particle *)p2)->r ) 
     return -1;
   
   return 0;
 }
void load_data(char *datafile)
{
    int i,j;
    FloatMat A;
    GenericMat B;
       
    printf("loading %s...\n",datafile);
    sprintf(A.name,"/x");
    load_hdfmatrixF(datafile,&A,1);
    if(A.size[1]!=3)
    {
      printf("Error, unexpected matrix size %d,%d\n",(int)A.size[0],(int)A.size[1]);
      exit(1);
    }
    nPall=A.size[0];
    Pall=malloc(sizeof(Particle)*nPall);
    for(i=0;i<nPall;i++)
    {
      for(j=0;j<3;j++)
      Pall[i].x[j]=A.x[i*3+j];
    }
    free(A.x);
    
    sprintf(A.name,"/v");
    load_hdfmatrixF(datafile,&A,1);
    if(A.size[1]!=3||A.size[0]!=nPall)
    {
      printf("Error, unexpected matrix size %zd,%zd (expecting %d,3)\n",A.size[0],A.size[1],nPall);
      exit(1);
    }
    for(i=0;i<nPall;i++)
    {
      for(j=0;j<3;j++)
      Pall[i].v[j]=A.x[i*3+j];
    }
    free(A.x);
    
    
    sprintf(B.name,"/flag");
    size_t nload;
    nload=load_hdfmatrix(datafile,&B,1,H5T_NATIVE_INT);
    if(nload==0) printf("Assuming flag=1 for every particle\n");
    int *p=B.x;
    for(i=0;i<nPall;i++)
      if(nload>0)
	Pall[i].flag=p[i];
      else
	Pall[i].flag=1;
      
      /*    
    double x0[3]={0.},v0[3]={0.};
    for(i=0;i<nPall;i++)
    {
      for(j=0;j<3;j++)
      {
	x0[j]+=Pall[i].x[j];
	v0[j]+=Pall[i].v[j];
      }
    }
    for(j=0;j<3;j++)
    {
      x0[j]/=nPall;
      v0[j]/=nPall;
    }
    for(i=0;i<nPall;i++)//shift to center
    {
      for(j=0;j<3;j++)
      {
	Pall[i].x[j]-=x0[j];
	Pall[i].v[j]-=v0[j];
      }
    }
*/    
    #define VecNorm(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    #define VecProd(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
    for(i=0;i<nPall;i++)
    {
      Pall[i].r=sqrt(VecNorm(Pall[i].x));
      Pall[i].K=VecNorm(Pall[i].v)/2.;
      Pall[i].L2=cross_product_norm2(Pall[i].x,Pall[i].v);
      Pall[i].vr=VecProd(Pall[i].x,Pall[i].v)/Pall[i].r; //radial vel
    }
    //sort
  // qsort(Pall, nPall, sizeof(Particle), cmpPartR); //sort according to r
    printf("%d particles loaded\n",nPall);
}

void sample_data(int subsample_id)
{//to be done: bootstrap when subsample_id<0
#ifdef SUBSAMPLE_SIZE
  nP=SUBSAMPLE_SIZE;
#else
  nP=nPall;
#endif
  P=malloc(sizeof(Particle)*nP);
  
  int i,j,imin,imax;
  imin=subsample_id*SUBSAMPLE_SIZE;
  imax=imin+SUBSAMPLE_SIZE;
  if(imin>=nPall) {printf("error: subsample id overflow; sampleid=%d,imin=%d,nP=%d\n", subsample_id, imin, nPall); exit(1);}
  if(imax>nPall) {printf("Warning: subsample truncated to tail of data; sampleid=%d, imax=%d, nP=%d\n", subsample_id, imax, nPall); imax=nPall;}
  for(i=imin,j=0;i<imax;i++)
  {
    if(Pall[i].r>R_MIN&&Pall[i].r<R_MAX
      #if defined(FILTER_FLAG)|defined(FILTER_RAND_SIZE)
      &&Pall[i].flag
      #endif
    )
    {
      P[j]=Pall[i];
      j++;
    }
  }
  nP=j;
  P=realloc(P,sizeof(Particle)*nP);
//   srand48(nP);
//   #ifdef FILTER_RAND_SIZE
//   for(i=0;i<nP;i++)
//     P[i].flag=(drand48()*nP<FILTER_RAND_SIZE);
//   #endif
}
