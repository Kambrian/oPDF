typedef struct
{
  int flag; //will be assigned a binID for binned MinDist estimator
  double r;
  double K;
  double L2;
  double x[3];
  double v[3];
  double E; //-(K+psi), binding energy; differ from previous version
  double T; //period=2*T
  double vr; //radial vel
  double theta; //radial phase, t/T
  double rlim[2]; //peri and apo center distance (estimated)
} Particle_t;
  
typedef struct 
{
  int nP;
  int nbin_r;
  int *RadialCount;
  double rmin, rmax;
  Particle_t *P;
} Tracer_t;

extern void load_tracer(char *datafile, Tracer_t * Sample);
extern void shuffle_tracer(unsigned long int seed, Tracer_t *Sample);
extern void resample_tracer(unsigned long int seed, Tracer_t *ReSample, Tracer_t *Sample);
extern void copy_tracer(int offset, int sample_size, Tracer_t *Sample, Tracer_t *FullSample);
extern void squeeze_tracer(Tracer_t *Sample);
extern void free_tracer(Tracer_t *Sample);
extern void print_tracer(Tracer_t *Sample);
extern void sort_tracer_flag(Tracer_t *Sample);
extern void cut_tracer(Tracer_t *Sample, double rmin, double rmax);
extern void count_tracer_radial(Tracer_t *Sample, int nbin);

extern void init_tracer(Tracer_t *Sample);
extern void make_sample(int sample_id, Tracer_t *Sample, Tracer_t *FullSample);