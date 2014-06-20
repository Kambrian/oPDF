#include <math.h>
#include "nr.h"
#define FUNC(x) (2.0*(x)*(*(funk->function))(bb-(x)*(x), funk->params))

double midsqu_gsl(gsl_function *funk, double aa, double bb, int n)
{
	double x,tnm,sum,del,ddel,a,b;
	static double s;
	#pragma omp threadprivate(s)
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC
