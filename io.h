#define R_MIN 0.
#define R_MAX 1000

#define Rhos0 2.187762e-3
#define Rs0 15.2
#define M0 200.
#define C0 15.
#define Z0 0.

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

extern void load_data(char *datafile);
