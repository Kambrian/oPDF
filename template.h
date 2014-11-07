#ifndef TEMPLATE_HEADER_INCLUDED

extern void init_potential_spline(int TMPid);
extern int get_current_TMPid();
extern void free_potential_spline();
extern double eval_potential_spline(double r);
extern double eval_density_spline(double r);
extern void decode_TemplateProf(Halo_t *halo);

#define TEMPLATE_HEADER_INCLUDED
#endif