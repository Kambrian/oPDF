#define R_MIN 40
#define R_MAX 50

typedef struct
  {
    int flag;
    double r;
    double K;
    double L2;
    double x[3];
    double v[3];
    double E; //K+psi
    double T; //period
    double rlim[2]; //peri and apo center distance (estimated)
  } Particle;
  
extern Particle *P;
extern int nP;

extern void load_data(char *datafile);
