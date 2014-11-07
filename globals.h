#ifndef GLOBALS_HEADER_INCLUDED

#define EPS 1e-16

/* virial definitions */
typedef enum
{ VIR_TH=0,
  VIR_C200,
  VIR_B200
} VirType_t;

struct global
{
  struct
  {
	double bin, bin_abs, rel;
  } tol;
  struct
  {
	double OmegaM0, OmegaL0;
  } cosmology;
  struct
  {
	double MassInMsunh, LengthInKpch, VelInKms;
	struct
	{
	  double G, H0;
	} Const;
  } units;
//   struct
//   {
// 	int views_NBinR;
// 	int views_RBinLog;
//   } misc;
};
extern struct global Globals;
extern void default_global_pars();
extern void set_units(double MassInMsunh, double LengthInKpch, double VelInKms);

#define GLOBALS_HEADER_INCLUDED
#endif