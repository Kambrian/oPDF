#ifndef MYMATH_HEADER_INCLUDED

#include <gsl/gsl_integration.h>
#ifndef INFINITY
  #define INFINITY (1.0/0.0)
#endif
#ifndef NAN 
  #define NAN (0.0/0.0)
#endif

#define POLY2(p1,p2,p3,x) ((p1*x+p2)*x+p3)
#define POLY5(p1,p2,p3,p4,p5,p6,x) (((((p1*x+p2)*x+p3)*x+p4)*x+p5)*x+p6)

#define MAXPAIR 64 //maximum number of pairs for a src gal.

typedef char Bool;
#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ fprintf(stderr,"Error opening file '%s'\n",filename);	fflush(stderr); exit(1);	}
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define VecNorm(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#define VecProd(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
#define SATURATE(x,x0) ((x)>(x0)?1.0:(x)/(x0))
#define BrokenPL(x,a,b) ((x)<1.0?pow(x,a):pow(x,b))
// debugging macros so we can pin down message origin at a glance
#define WHERESTR  "[file %s, line %d]: "
#define WHEREARG  __FILE__, __LINE__
#define DEBUGPRINT2(...)       fprintf(stderr, __VA_ARGS__)
#define DEBUGPRINT(_fmt, ...)  DEBUGPRINT2(WHERESTR _fmt, WHEREARG, __VA_ARGS__)
//example:  DEBUGPRINT("hey, x=%d\n", x);


extern int cmpDouble(const void *p1, const void *p2);
extern double factorial(int n);
extern double log_factorial(int n);
extern double smpintD(gsl_function *F, double lowerlim,double upperlim, double tol);
extern int Dmax_of_vec(double *vec,int Len);//what if the max is not unique???????
extern int Fmax_of_vec(float *vec,int Len);//what if the max is not unique?,return first max
extern int Dmin_of_vec(double *vec,int Len);//what if the max is not unique???????
extern int Fmin_of_vec(float *vec,int Len);//what if the max is not unique?,return first min
extern int insert_to_array(int a, int len, int *ascd_arr);//insert into an ascending arr
extern int check_dup(int l1, int *p1,int l2, int *p2);  //check duplicating IDs for two arrays
extern float distance(float x[3],float y[3]);
extern double distanceD(double x[3],double y[3]);
extern double cross_product_norm2(double a[3],double b[3]);
extern double cosdistD(double x[2], double y[2]);
extern void *mymalloc(size_t n);
extern void myfree(void *mem);
extern void spherical_basisD(double dx[3],double er[3],double et[3],double ef[3]);//er: radial; et: azimuthal,(0~2*pi); ef: elevation (0~pi) 
extern float f_prod(float *a,float *b,int dim);
extern void f_cross(float a[3],float b[3],float c[3]);

extern size_t memtrim(void *mem,size_t membersize,char *mask,size_t n);
extern size_t memselect(void * memto, const void * memfrom, size_t membersize, char *mask, size_t n);
extern void * mem_icopy(const void * memfrom, int *idlist, int nid, size_t membersize);
extern void init_progress_monitor(int ini, int totstep, int nprint);
// ini: initial value progress percentage, set to 0 if no reason
//totstep: maximum iteration steps to monitor
//nprint: number of revealing times
extern void monitor_progress(int istep);
//put this inside loop to monitor progress.
extern int try_readfile(char * filename);
extern double KSprob(int n,double d);
extern float KSprobNR(int n, double d);
extern float KuiperProb(int n, double d);
extern double GeneralizedExtremeValuePDF(double x, double mu, double sigma, double k);
extern double NormPDF(double x,double mu, double sigma);
#define MYMATH_HEADER_INCLUDED
#endif
