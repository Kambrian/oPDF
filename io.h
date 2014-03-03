typedef struct
  {
    double r;
    double K;
    double L2;
    double x[3];
    double v[3];
  } Particle;
  
extern Particle *P;
extern int nP;

extern void load_data(char *datafile);
