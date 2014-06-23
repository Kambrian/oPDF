// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/
#include <iostream>
#include "Minuit2/FCNGradientBase.h"
extern "C" {
#include "io.h"
#include "models.h"
}

namespace ROOT {

   namespace Minuit2 {


class TracerLike : public FCNBase {

public:
  Tracer_t *Sample;
  int Estimator;
  
  TracerLike(Tracer_t& sample, int estimator): Sample(&sample), Estimator(estimator), errdef(1.0) {
  }

  ~TracerLike() {}

  double operator()(const std::vector<double>& par) const {
// 	std::cout<<par[0]<<","<<par[1]<<std::endl;
    return -freeze_and_like(&par[0], Estimator, Sample);
  }
  
  void SetErrorDef(double err) { errdef=err;}
  double Up() const {return errdef;}

private:
  double errdef;
};

  }  // namespace Minuit2

}  // namespace ROOT
