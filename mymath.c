#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "mymath.h"

int cmpDouble(const void *p1, const void *p2)
{ //in ascending order
  if((*(double *)p1) > (*(double *)p2) ) 
    return 1;
  
  if((*(double *)p1) < (*(double *)p2) ) 
    return -1;
  
  return 0;
}
double factorial(int n)
{
  double result=1.;

  while(n>0)
  {
    result*=n;
    n--;
  }

  return result;
}
double log_factorial(int n)
{
  double result=0.;

  while(n>0)
  {
    result+=log(n);
    n--;
  }

  return result;
}
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

//KS test code from http://www.jstatsoft.org/v08/i18/paper
static void mMultiply(double *A,double *B,double *C,int m)
{int i,j,k; double s;
for(i=0;i<m;i++) for(j=0; j<m; j++)
{s=0.; for(k=0;k<m;k++) s+=A[i*m+k]*B[k*m+j]; C[i*m+j]=s;}
}
static void mPower(double *A,int eA,double *V,int *eV,int m,int n)
{ double *B;int eB,i;
if(n==1) {for(i=0;i<m*m;i++) V[i]=A[i];*eV=eA; return;}
mPower(A,eA,V,eV,m,n/2);
B=(double*)malloc((m*m)*sizeof(double));
mMultiply(V,V,B,m);
eB=2*(*eV);
if(n%2==0){for(i=0;i<m*m;i++) V[i]=B[i]; *eV=eB;}
else {mMultiply(A,B,V,m);
*eV=eA+eB;}
if(V[(m/2)*m+(m/2)]>1e140) {for(i=0;i<m*m;i++) V[i]=V[i]*1e-140;*eV+=140;}
free(B);
}
double KSprob(int n,double d)
{ //return K(n, d) = Pr(Dn >= d), probability of compatibility
  int k,m,i,j,g,eH,eQ;
double h,s,*H,*Q;
//OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 2.*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
k=(int)(n*d)+1; m=2*k-1; h=k-n*d;
H=(double*)malloc((m*m)*sizeof(double));
Q=(double*)malloc((m*m)*sizeof(double));
for(i=0;i<m;i++) for(j=0;j<m;j++)
if(i-j+1<0) H[i*m+j]=0; else H[i*m+j]=1;
for(i=0;i<m;i++) {H[i*m]-=pow(h,i+1); H[(m-1)*m+i]-=pow(h,(m-i));}
H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
for(i=0;i<m;i++) for(j=0;j<m;j++)
if(i-j+1>0) for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
eH=0; mPower(H,eH,Q,&eQ,m,n);
s=Q[(k-1)*m+k-1];
for(i=1;i<=n;i++) {s=s*i/n; if(s<1e-140) {s*=1e140; eQ-=140;}}
s*=pow(10.,eQ); free(H); free(Q); return 1.-s;
}

#define EPS_ITER 1e-3
#define EPS_TOT 1.0e-8
float KSprobNR(int n, double d)
{//from numerical recipes, P(D>Dobs), tail probability
  float en=sqrt(n);
  float lambda2=d*(en+0.12+0.11/en);
  lambda2*=-2.*lambda2;
  int j;
  float fac=2.0,sum=0.0,term,termbf=0.0;
  
  for (j=1;j<=100;j++) 
  {
    term=fac*exp(lambda2*j*j);
    sum += term;
    term=fabs(term);
    if (term <= EPS_ITER*termbf || term <= EPS_TOT*sum) return sum;
    fac = -fac;// Alternating signs in sum.
    termbf=term;
  }
  return 1.0; //Get here only by failing to converge.
}

float KuiperProb(int n, double d)
{//ref: numerical recipes
  float en=sqrt(n);
  float lambda2=(en+0.155+0.24/en)*d;
  if(lambda2<0.4) return 1.0;
  
  lambda2*=lambda2*2;
  int j;
  float xj,sum=0.0,term,termbf=0.0;
  
  for (j=1;j<=100;j++) 
  {
    xj=lambda2*j*j;
    term=2.*(2.*xj-1.)*exp(-xj);
    sum += term;
    term=fabs(term);
    if (term <= EPS_ITER*termbf || term <= EPS_TOT*sum) return sum;
    termbf=term;
  }
  return 1.0;
}

double GeneralizedExtremeValuePDF(double x, double mu, double sigma, double k)
{
  double t, tk;
  if(k==0.)
  {
	t=-(x-mu)/sigma;
	return 1./sigma*exp(t-exp(t));
  }
  
  t=1+k*(x-mu)/sigma;
  if(t<=0) return 0.;
  tk=pow(t, -1./k);
  return 1./sigma*exp(-tk)*tk/t;
}

double NormPDF(double x,double mu, double sigma)
{
	return exp(-(x-mu)*(x-mu)/sigma/sigma/2.)/sqrt(2*M_PI)/sigma;
}