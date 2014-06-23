//playing around with machineprecision is still no good; dangerous. stick to scipy.fmin() please.
//!!IMPORTANT: --disable-openmp when compile MINUIT2, to avoid conflicting with user-func's openmp here.
#include "TracerMinuit.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"

// StackAllocator gStackAllocator;

using namespace ROOT::Minuit2;
// using namespace std;

int main() {

  int samplesize=1000, estimator=8;
  Tracer_t FullSample={}, Sample={}; //initialize to zeros!!!
  init_tracer(&FullSample);
  make_sample(0, samplesize, &Sample, &FullSample);
  free_tracer(&FullSample);
  
  alloc_integration_space();
  
  TracerLike fcn(Sample, estimator);
  fcn.SetErrorDef(1);
  
  MnUserParameters upar;
  upar.Add("m", 2., 0.1);
  upar.Add("c", 1., 0.1);
//   upar.SetPrecision(1e-4); //either set here for parameter-specific, or with minimizer.SetPrecision() for minimizer-specific precision
  std::cout<<upar.Precision()<<std::endl;
  std::cout<<upar<<std::endl;
  
//   MODEL_TOL_REL=1e-5;
  MnMigrad migrad(fcn, upar);
  migrad.SetPrecision(1e-4);
  FunctionMinimum min = migrad(200, 10.);
  std::cout<<"migrad minimum: "<<min<<min.UserCovariance()<<std::endl;  
  std::cout<<"--------------------------------"<<std::endl;
  
//   MnContours contours(fcn, min);
//   fcn.SetErrorDef(1);
//   std::vector<std::pair<double,double> > cont = contours(0, 1, 20);
//   MnPlot plot;
//   plot(min.UserState().Value("m"),min.UserState().Value("c"),cont);

//   MnHesse hesse; 
//   hesse( fcn, min); 
//   std::cout<<"minimum after hesse: "<<min<<std::endl;
  
//   MODEL_TOL_REL=1e-7;
  MnSimplex simplex(fcn, upar);
//   simplex.SetPrecision(1e-4);
  std::cout<<simplex.Precision()<<std::endl;
  FunctionMinimum smin = simplex(100, 1e-3);
  std::cout<<"simplex minimum: "<<smin<<std::endl;
  
  free_tracer(&Sample);
  free_integration_space();
  
  return 0;
}
