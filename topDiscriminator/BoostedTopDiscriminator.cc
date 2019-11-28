#include <cmath>
#include "BoostedTopDiscriminator.h"

BoostedTopDiscriminator::BoostedTopDiscriminator(TString weights)
{
    this->weights = weights;
    reader_ = new TMVA::Reader("!Color:!Silent");

    reader_->AddVariable("jetTau3", &var_[0]);
    reader_->AddVariable("jetTau2", &var_[1]);
    reader_->AddVariable("jetTau1", &var_[2]);
    reader_->AddVariable("jetMassSub0", &var_[3]);
    reader_->AddVariable("jetMassSub1", &var_[4]);

    reader_->BookMVA("MVA", weights);
}

float BoostedTopDiscriminator::eval(float tau3, float tau2, float tau1, float jetMassSub0, float jetMassSub1)
{
  var_[0] = tau3;
  var_[1] = tau2;
  var_[2] = tau1;
  var_[3] = jetMassSub0;
  var_[4] = jetMassSub1;

  return reader_->EvaluateMVA("MVA");
}

BoostedTopDiscriminator::~BoostedTopDiscriminator()
{
    delete reader_;
}
