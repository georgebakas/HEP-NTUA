#ifndef BoostedTopCategoryDiscriminator1_h
#define BoostedTopCategoryDiscriminator1_h

#include "TString.h"
#include "TMVA/Reader.h"

class BoostedTopCategoryDiscriminator1
{
  public:
    BoostedTopCategoryDiscriminator1(TString weights);
    ~BoostedTopCategoryDiscriminator1();
    float eval(float tau3, float tau2, float tau1, float jetMassSub0, float jetMassSub1, 
    			float ecfB1N2, float ecfB1N3, float ecfB2N2, float ecfB2N3, float JetPtOverSumPt, float jetPt);

  private:
    TString weights;
    TMVA::Reader *reader_;
    float var_[11];
};

#endif
