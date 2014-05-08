#define DATAFILE "mockhalo_wenting.hdf5"

#define R_MIN 1
#define R_MAX 1000.
// #define FILTER_FLAG
// #define FILTER_RAND_SIZE 10000 //select a random sample of the given size

#define SUBSAMPLE_SIZE 1000 //only load a subsample of this size 

/* mockhalo--*/
#define Rhos0 2.187762e-3
#define Rs0 15.2
#define M0  183.5017
#define C0 16.1560
// */
#define Z0 0.

/*A4--
#define Rhos0 0.00220633
#define Rs0 15.1576
#define M0 183.8
#define C0 16.21
// */
/*B4--*
#define Rhos0 0.000515264
#define Rs0 20.9363
#define M0 83.45
#define C0 9.02
// */

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
  
extern Particle *P;
extern int nP;

extern void load_data(char *datafile, int dataid);
extern void sample_data(int subsample_id);