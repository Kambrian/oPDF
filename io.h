// #define FILTER_FLAG
// #define FILTER_RAND_SIZE 10000 //select a random sample of the given size

// #define SUBSAMPLE_SIZE 500000 //only load a subsample of this size 

// #define V0
// #define L0
// #define E0

typedef struct
  {
    int flag;
    double r;
    double K;
    double L2;
    double x[3];
    double v[3];
    double E; //K+psi
    double T; //period=2*T
    double vr; //radial vel
    double theta; //radial phase, t/T
    double rlim[2]; //peri and apo center distance (estimated)
  } Particle;

extern double R_MIN, R_MAX;  
extern Particle *P;
extern int nP;

extern void load_data(char *datafile);
extern void shuffle_data(unsigned long int seed);
extern void sample_data(int subsample_id);
extern int squeeze_data();
extern void free_sample();
extern void free_data();