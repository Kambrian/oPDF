#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mymath.h"
#include "hdf_util.h"
#include "io.h"
  
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
    sprintf(A.name,"/PartType0/Coordinates");
    load_hdfmatrixF(datafile,&A,1);
    if(A.size[1]!=3)
    {
      printf("Error, unexpected matrix size %d,%d\n",(int)A.size[0],(int)A.size[1]);
      exit(1);
    }
    nP=A.size[0];
    P=malloc(sizeof(Particle)*nP);
    for(i=0;i<nP;i++)
    {
      for(j=0;j<3;j++)
      P[i].x[j]=A.x[i*3+j];
    }
    free(A.x);
    
    sprintf(A.name,"/PartType0/Velocities");
    load_hdfmatrixF(datafile,&A,1);
    if(A.size[1]!=3||A.size[0]!=nP)
    {
      printf("Error, unexpected matrix size %zd,%zd (expecting %d,3)\n",A.size[0],A.size[1],nP);
      exit(1);
    }
    for(i=0;i<nP;i++)
    {
      for(j=0;j<3;j++)
      P[i].v[j]=A.x[i*3+j];
    }
    free(A.x);
    
    sprintf(B.name,"/PartType0/Flag");
    load_hdfmatrix(datafile,&B,1,H5T_NATIVE_INT);
    int *p=B.x;
    for(i=0;i<nP;i++)
      P[i].flag=p[i];
    
//     printf("%g,%g,%g;%g,%g,%g\n", P[5].x[0],P[5].x[1],P[5].x[2], P[10].v[0],P[10].v[1],P[10].v[2]);
/*    
    double x0[3]={0.},v0[3]={0.};
    for(i=0;i<nP;i++)
    {
      for(j=0;j<3;j++)
      {
	x0[j]+=P[i].x[j];
	v0[j]+=P[i].v[j];
      }
    }
    for(j=0;j<3;j++)
    {
      x0[j]/=nP;
      v0[j]/=nP;
    }
    for(i=0;i<nP;i++)//shift to center
    {
      for(j=0;j<3;j++)
      {
	P[i].x[j]-=x0[j];
	P[i].v[j]-=v0[j];
      }
    }
*/    
    #define VecNorm(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    #define VecProd(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
    for(i=0;i<nP;i++)
    {
      P[i].r=sqrt(VecNorm(P[i].x));
      P[i].K=VecNorm(P[i].v)/2.;
      P[i].L2=cross_product_norm2(P[i].x,P[i].v);
      P[i].vr=VecProd(P[i].x,P[i].v)/P[i].r; //radial vel
    }
//     printf("%g,%g,%g\n",P[2].r,P[2].K,P[2].L2);
    //select particles....
    for(i=0,j=0;i<nP;i++)
    {
      if(P[i].r>R_MIN&&P[i].r<R_MAX&&P[i].flag)
      {
	if(i>j)	memmove(P+j,P+i,sizeof(Particle));
	j++;
      }
    }
    nP=j;
    P=realloc(P,sizeof(Particle)*nP);
    //sort
   // qsort(P, nP, sizeof(Particle), cmpPartR); //sort according to r
    printf("%d particles loaded\n",nP);
} 
