#include <cmath>
#include "BoostedTopCategoryDiscriminator_tTagger.h"

BoostedTopCategoryDiscriminator_tTagger::BoostedTopCategoryDiscriminator_tTagger(TString weights)
{
    this->weights = weights;
    reader_ = new TMVA::Reader("!Color:!Silent");

    reader_->AddVariable("jetTau3", &var_[0]);
    reader_->AddVariable("jetTau2", &var_[1]);
    reader_->AddVariable("jetTau1", &var_[2]);
    reader_->AddVariable("jetMassSub0", &var_[3]);
    reader_->AddVariable("jetMassSub1", &var_[4]);

    reader_->AddVariable("ecfB1N2", &var_[5]);
    reader_->AddVariable("ecfB1N3", &var_[6]);
    reader_->AddVariable("ecfB2N2", &var_[7]);
    reader_->AddVariable("ecfB2N3", &var_[8]);
    reader_->AddVariable("JetPtOverSumPt", &var_[9]);

    reader_->AddSpectator("jetPt", &var_[10]);

    reader_->BookMVA("BDTCat", weights);
}

float BoostedTopCategoryDiscriminator_tTagger::eval(float tau3, float tau2, float tau1, float jetMassSub0, float jetMassSub1, 
                                             float ecfB1N2, float ecfB1N3, float ecfB2N2, float ecfB2N3, float JetPtOverSumPt, 
                                             float jetPt)
{
  var_[0] = tau3;
  var_[1] = tau2;
  var_[2] = tau1;
  var_[3] = jetMassSub0;
  var_[4] = jetMassSub1;
  var_[5] = ecfB1N2;
  var_[6] = ecfB1N3;
  var_[7] = ecfB2N2;
  var_[8] = ecfB2N3;
  var_[9] = JetPtOverSumPt;
  var_[10] = jetPt;

  return reader_->EvaluateMVA("BDTCat");
}

BoostedTopCategoryDiscriminator_tTagger::~BoostedTopCategoryDiscriminator_tTagger()
{
    delete reader_;
}
