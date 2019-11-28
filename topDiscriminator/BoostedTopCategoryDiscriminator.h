#ifndef BoostedTopCategoryDiscriminator_h
#define BoostedTopCategoryDiscriminator_h

#include "TString.h"
#include "TMVA/Reader.h"

class BoostedTopCategoryDiscriminator
{
  public:
    BoostedTopCategoryDiscriminator(TString weights);
    ~BoostedTopCategoryDiscriminator();
    float eval(float tau3, float tau2, float tau1, float jetMassSub0, float jetMassSub1, float jetPt);

  private:
    TString weights;
    TMVA::Reader *reader_;
    float var_[6];
};

#endif
