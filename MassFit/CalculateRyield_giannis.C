float CalculateRyield(TString variable, float &error)
{
  Double_t rYield, rYield_error;
  //reduced
  Double_t mcSig_2btag_reduced_error, mcSig_0btag_reduced_error, data_0btag_reduced_error;
  Double_t mcSig_2btag_reduced = GetHistoIntegralAndError(AnalysisConstants::files[year]["mcBkg"], TString::Format("2btag_%s", variable.Data()), AnalysisConstants::luminositiesSR[year], mcSig_2btag_reduced_error);
  //std::cout << "mcSig 0btag reduced: ";
  Double_t mcSig_0btag_reduced = GetHistoIntegralAndError(AnalysisConstants::files[year]["mcBkg"], TString::Format("0btag_%s", variable.Data()), AnalysisConstants::luminositiesCR[year], mcSig_0btag_reduced_error);
  //std::cout << "data 0btag reduced: ";
  Double_t data_0btag_reduced = GetHistoIntegralAndError(AnalysisConstants::files[year]["data"], TString::Format("0btag_%s", variable.Data()), AnalysisConstants::luminositiesCR[year], data_0btag_reduced_error);

  //extended
  Double_t mcSig_2btag_extended_error, mcSig_0btag_extended_error, data_0btag_extended_error;
  Double_t mcSig_2btag_extended = GetHistoIntegralAndError(AnalysisConstants::files[year]["mcBkgExtended"], TString::Format("2btag_%s", variable.Data()), AnalysisConstants::luminositiesSR[year], mcSig_2btag_extended_error);
  //std::cout << "mcSig 0btag extended: ";
  Double_t mcSig_0btag_extended = GetHistoIntegralAndError(AnalysisConstants::files[year]["mcBkgExtended"], TString::Format("0btag_%s", variable.Data()), AnalysisConstants::luminositiesCR[year], mcSig_0btag_extended_error);
  //std::cout << "data 0btag extended: ";
  Double_t data_0btag_extended = GetHistoIntegralAndError(AnalysisConstants::files[year]["dataExtended"], TString::Format("0btag_%s", variable.Data()), AnalysisConstants::luminositiesCR[year], data_0btag_extended_error);

  rYield = (mcSig_2btag_reduced / mcSig_2btag_extended) / (mcSig_0btag_reduced / mcSig_0btag_extended);
  //std::cout << "rYield initial: " << rYield << std::endl;
  Double_t rYield_correction = (data_0btag_reduced / data_0btag_extended);
  rYield = rYield * rYield_correction;
  //std::cout << "rYield corrected: " << rYield << std::endl;

  rYield_error = TMath::Sqrt(TMath::Power(mcSig_0btag_extended / (mcSig_0btag_reduced * mcSig_2btag_extended)  mcSig_2btag_reduced_error, 2) +
                             TMath::Power(mcSig_2btag_reduced / (mcSig_0btag_reduced  * mcSig_2btag_extended)  mcSig_0btag_extended_error, 2) +
                             TMath::Power((mcSig_2btag_reduced * mcSig_0btag_extended) / (TMath::Power(mcSig_0btag_reduced, 2) * mcSig_2btag_extended) * mcSig_0btag_reduced_error, 2) +
                             TMath::Power((mcSig_2btag_reduced * mcSig_0btag_extended) / (mcSig_0btag_reduced * TMath::Power(mcSig_2btag_extended, 2)) * mcSig_2btag_extended_error, 2));

  Double_t rYield_correction_error = TMath::Sqrt(TMath::Power(1 / data_0btag_extended * data_0btag_reduced_error, 2) +
                                                 TMath::Power(data_0btag_reduced / TMath::Power(data_0btag_extended, 2) * data_0btag_extended_error, 2));

  error = (float)TMath::Sqrt(TMath::Power(rYield_correction * rYield_error, 2) + TMath::Power(rYield * rYield_correction_error, 2));
  std::cout << "Error: " << error << std::endl;
  return (float)rYield;
}
