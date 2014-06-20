#include <gsl/gsl_math.h>
#ifndef _NR_H_
#define _NR_H_
extern void nrerror(char error_text[]);
extern double qromo_gsl(gsl_function *func, double a, double b,
			 double (*choose)(gsl_function *, double, double, int), double tol_rel);
extern double midsql_gsl(gsl_function *funk, double aa, double bb, int n);
extern double midsqu_gsl(gsl_function *funk, double aa, double bb, int n);

#endif /* _NR_H_ */
