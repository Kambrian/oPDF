struct SplineData
{
	int FlagUseSpline; //spline ready to be used. override the default potential calculation with spline interpolation
	gsl_interp_accel *acc;
	gsl_spline *spline ;
	gsl_interp_accel *acc_dens;
	gsl_spline *spline_dens;
	int TMPid;
	double Rs;
};

extern void init_potential_spline(int TMPid);
extern void free_potential_spline();
extern double eval_potential_spline(double r);
extern double eval_density_spline(double r);
extern void decode_TemplateProf(Halo_t *halo);
