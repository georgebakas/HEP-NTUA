#include <cmath>
#include "BoostedTopCategoryDiscriminator.h"

BoostedTopCategoryDiscriminator::BoostedTopCategoryDiscriminator(TString weights)
{
    this->weights = weights;
    reader_ = new TMVA::Reader("!Color:!Silent");

    reader_->AddVariable("jetTau3", &var_[0]);
    reader_->AddVariable("jetTau2", &var_[1]);
    reader_->AddVariable("jetTau1", &var_[2]);
    reader_->AddVariable("jetMassSub0", &var_[3]);
    reader_->AddVariable("jetMassSub1", &var_[4]);
    
    reader_->AddSpectator("jetPt", &var_[5]);

    reader_->BookMVA("BDTCat", weights);
}

float BoostedTopCategoryDiscriminator::eval(float tau3, float tau2, float tau1, float jetMassSub0, float jetMassSub1, float jetPt)
{
  var_[0] = tau3;
  var_[1] = tau2;
  var_[2] = tau1;
  var_[3] = jetMassSub0;
  var_[4] = jetMassSub1;
  //std::cout<<"Pt: "<<jetPt<<std::endl;
  var_[5] = jetPt;

  return reader_->EvaluateMVA("BDTCat");
}

BoostedTopCategoryDiscriminator::~BoostedTopCategoryDiscriminator()
{
    delete reader_;
}
