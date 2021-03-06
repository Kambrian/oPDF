#ifndef TRACER_HEADER_INCLUDED

#include "halo.h"

typedef struct
{
  int haloid;//host haloid
  int subid;//subhalo id
//   int strmid; //stream id; 11/09/2014:suppressed.
  int flag; //will be assigned a binID for binned MinDist estimator
  double w; //weight (particle mass normalized by average mass), mostly for tagged stars.
  double r;
  double K;
  double L2;//L^2
  double L;//L
  double x[3];
  double v[3];
  double E; //-(K+psi), binding energy; differ from previous version
  double T; //period=2*T
  double vr; //radial vel
  double theta; //radial phase, t/T
  double rlim[2]; //peri and apo center distance (estimated)
} Particle_t;

struct Tracer;
typedef struct Tracer Tracer_t;
typedef Tracer_t TracerView; //the P pointer are not allocated, but just points to other Tracer's data; so they are only views.

struct Tracer
{
  double lnL; //the likelihood associated with the lowest level view through like_eval(); upper level views are chi2.
  int nP;
  double mP; //the average mass of particles
  Particle_t *P;
  int nbin_r; //never set this manually. only modify it via count_tracer_radial()
  int FlagRLogBin; //whether to bin r in logspace or linear space
  double *RadialCount;
  double rmin, rmax;
  double proxybin[2]; //for the edges of views.
  Halo_t halo;
  int nView;
  char ViewType;
  TracerView *Views;
};

extern void load_tracer_particles(char *datafile, char *grpname, Tracer_t * Sample, int AddHubbleFlow);
extern void cut_tracer_particles(Tracer_t *Sample, double rmin, double rmax);
extern void shuffle_tracer_particles(unsigned long int seed, Tracer_t *Sample);
extern void squeeze_tracer_particles(Tracer_t *Sample);
extern void print_tracer_particle(Tracer_t *Sample, int i);
extern void free_tracer_particles(Tracer_t *Sample);

extern void copy_tracer_particles(int offset, int sample_size, Tracer_t *Sample, Tracer_t *FullSample);
extern void resample_tracer_particles(unsigned long int seed, Tracer_t *ReSample, Tracer_t *Sample);

extern void count_tracer_radial(Tracer_t *Sample, int nbin, int FlagRLogBin);
extern void free_tracer_rcounts(Tracer_t * Views);

extern void sort_part_flag(Particle_t *P, int nP);
extern void sort_part_L(Particle_t *P, int nP);
extern void sort_part_E(Particle_t *P, int nP);
extern void create_tracer_views(Tracer_t *Sample, int nView, char proxy);
extern void free_tracer_views(TracerView * Views);
extern void create_nested_views(int nbin[], char ViewTypes[],  Tracer_t *Sample);

extern void free_tracer(Tracer_t *Sample);

void tracer_set_energy(Tracer_t *Sample);
void tracer_set_orbits(Tracer_t *Sample, int FlagSetPhase);

#define TRACER_HEADER_INCLUDED
#endif
