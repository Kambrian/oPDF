#include <math.h>
#include <stdio.h>
#include "nr.h"
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

double qromo_gsl(gsl_function *func, double a, double b,
	double (*choose)(gsl_function *, double, double, int), double tol_rel)
{
	int j;
	double ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < tol_rel*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
// 	fprintf(stderr, "Warning: Too many steps in routing qromo.");
	return ss;
}
#undef JMAX
#undef JMAXP
#undef K
