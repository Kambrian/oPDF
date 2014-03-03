#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "mymath.h"

double smpintD(gsl_function *F, double lowerlim,double upperlim, double tol)
{
	//double tol=1e-3; //relative tolerance
	int maxit=15,N0=8,i,count=1; //first divide into 2*N0 intervals
	
	double sum0,sum1,sum2,sum3,sint;
	double step2,x,xrng[2];
	
	xrng[0]=lowerlim;
	xrng[1]=upperlim;
	if(xrng[0]==xrng[1])
	return 0;
	
	#define rfun(x) F->function(x,F->params)
	
	sum0=rfun(xrng[0])+rfun(xrng[1]);//the two ends
	step2=(xrng[1]-xrng[0])/N0;					//double of the step
	sum1=sum2=0;
	for(i=0,x=xrng[0]+step2/2;i<N0;i++,x+=step2) //odd points
	{
		sum1+=rfun(x);		
	}
	for(i=0,x=xrng[0]+step2;i<N0-1;i++,x+=step2)	//even points
	{
		sum2+=rfun(x);		
	}
	
	//new odd points
	while(count<maxit)
	{
		sum3=0;N0*=2;step2/=2;
		for(i=0,x=xrng[0]+step2/2;i<N0;i++,x+=step2)
		{
			sum3+=rfun(x);
		}
		sint=(sum0+2*(sum1+sum2)+4*sum3);
		if(fabs(-sum0-6*sum1-2*sum2+4*sum3)<tol*sint)
		{
			//printf("iter %d, div=%d\n", count, N0);
			return sint*step2/6;
		}
		sum2=sum1+sum2;
		sum1=sum3;
		count++;
		//printf("iter %d, div=%d\n", count, N0);
	}

// 	printf("warning: maximum number of iterations %d exceeded for relative error %f,integral stops here\n",maxit,tol);
	return sint*step2/6;
}

	
float distance(float x[3],float y[3])
{
	float dx[3];
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	//~ #ifdef PERIODIC_BDR
	//~ dx[0]=NEAREST(dx[0]);
	//~ dx[1]=NEAREST(dx[1]);
	//~ dx[2]=NEAREST(dx[2]);
	//~ #endif
	return sqrtf(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}
double distanceD(double x[3],double y[3])
{
	double dx[3];
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}
size_t memtrim(void *mem,size_t membersize,char *mask,size_t n)
{//keep all the members in mem with TRUE mask and move them to the beginning of mem
	return memselect(mem,mem,membersize,mask,n);
}
size_t memselect(void * memto, const void * memfrom, size_t membersize, char *mask, size_t n)
{//copy elements from memfrom to memto accordding to flags in mask;
size_t i,j;
char *pto,*pfrom;
pto=(char *)memto;
pfrom=(char *)memfrom;
j=0;
for(i=0;i<n;i++)
{
	if(mask[i])
	{
		memmove(pto+j*membersize,pfrom+i*membersize,membersize);
		j++;
	}
}
return j;	
}

void * mem_icopy(const void * memfrom, int *idlist, int nid, size_t membersize)
{/*to copy members specified in idlist 
* from memfrom 
* to a newly allocated array (memto)
*  nid is length of idlist
*  membersize is membersize of memfrom and memto
* return memto
*/
int i;
char *pto,*pfrom;
pto=malloc(membersize*nid);
pfrom=(char *)memfrom;
for(i=0;i<nid;i++)
	memcpy(pto+i*membersize,pfrom+idlist[i]*membersize,membersize);
return (void *)pto;
}

void spherical_basisD(double dx[3],double er[3],double et[3],double ef[3])
/*er: radial;
 *et: azimuthal,(0~2*pi)
 *ef: elevation (0~pi)
 * */ 
{
	int i;
	double dr,dxy2,modet,modef;//unit vectors of spherical coord.
	dxy2=dx[0]*dx[0]+dx[1]*dx[1];
	dr=sqrt(dxy2+dx[2]*dx[2]);
	for(i=0;i<3;i++) er[i]=dx[i]/dr;
	modet=sqrt(dxy2);
	if(0.==modet)
	{
		et[0]=1.;
		et[1]=et[2]=0.;
		ef[1]=1.;
		ef[0]=ef[2]=0.;
	}
	else
	{
		et[0]=-dx[1]/modet;
		et[1]=dx[0]/modet;
		et[2]=0;
		if(0.==dx[2])
		{
		ef[0]=ef[1]=0.;
		ef[2]=1.;	
		}
		else
		{
		modef=sqrt(dxy2+dxy2*dxy2/dx[2]/dx[2]);
		ef[0]=-dx[0]/modef;
		ef[1]=-dx[1]/modef;
		ef[2]=dxy2/dx[2]/modef;
		}
	}
}

float f_prod(float *a,float *b,int dim)
{
	int i;
	float c;
	c=0.;
	for(i=0;i<dim;i++)
		c+=a[i]*b[i];
	return c;
}
void f_cross(float a[3],float b[3],float c[3])
{
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[0]*b[2]-a[2]*b[0];
	c[2]=a[0]*b[1]-a[1]*b[0];
}

double cross_product_norm2(double a[3],double b[3])
{//return c**2
  double c[3];
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[0]*b[2]-a[2]*b[0];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return (c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
}

static int progress=0;
static int totsteps=0;
static int nshow=10;
void init_progress_monitor(int ini, int totstep, int nprint)
{// ini: initial value progress percentage, set to 0 if no reason
//totstep: maximum iteration steps to monitor
//nprint: number of revealing times
	progress=ini;
	totsteps=totstep;
	nshow=nprint;
	printf(" %02d%%",progress*100/nshow);fflush(stdout);
}
void monitor_progress(int istep)
{  //put this inside loop to monitor progress.
	if(istep>totsteps*progress/nshow)
	{
		printf("\b\b\b%02d%%",progress*100/nshow);fflush(stdout);
		//~ printf(",%02d%%",progress*100/nshow);fflush(stdout);
		progress++;
	}
}

int try_readfile(char * filename)
{// return 1 when file can be read; 0 otherwise.
	FILE *fp;
	if((fp=fopen(filename,"r")))
	{
		 fclose(fp);
		 return 1;
	 }
	return 0;
}

int seek_new_line(FILE *fp)
{
	char tmp;
	tmp=fgetc(fp);
	while(tmp!='\n')
	{
		if(feof(fp))
			return -1;
		tmp=fgetc(fp);
	}
//	LINE_NUM++;
	return 1;
}
int skip_comment(FILE *fp,char c)
{
//skip those lines starting with the specified char
	while(c==(fgetc(fp)))
	{
		if(feof(fp))
			return -1;
		seek_new_line(fp);
	}
	fseek(fp,-sizeof(char),SEEK_CUR);
	return 1;
}