#ifndef BoostedTopDiscriminator_h
#define BoostedTopDiscriminator_h

#include "TString.h"
#include "TMVA/Reader.h"

class BoostedTopDiscriminator
{
  public:
    BoostedTopDiscriminator(TString weights);
    ~BoostedTopDiscriminator();
    float eval(float tau3, float tau2, float tau1, float jetMassSub0, float jetMassSub1);

  private:
    TString weights;
    TMVA::Reader *reader_;
    float var_[5];
};

#endif
