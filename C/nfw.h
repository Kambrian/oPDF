#ifndef NFW_HEADER_INCLUDED

#include "tracer.h"

extern double NFW_concentration(double M, double z, int virtype);
extern void decode_NFWprof(Halo_t *halo);//translate halo (M,c) to other profile parameters
extern double NFW_mass(double r, Halo_t *halo);
extern double NFW_like(double pars[],Tracer_t *T);

#define NFW_HEADER_INCLUDED
#endif