#ifndef ANALYSISCONSTANTS_UL_H_
#define ANALYSISCONSTANTS_UL_H_

#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TStyle.h>

#include <map>
#include <vector>
#include <utility>
#include <iostream>

namespace AnalysisConstants_UL
{
  //std::vector<TString> years = {"2018"};
  std::vector<TString> years;
  bool debug = true;

  std::map<TString, float> floatConstants;
  std::map<TString, float> luminositiesSR, luminositiesCR;

  std::map<TString, float> fitConstants;
  std::map<TString, float> fitConstantsErrors;

  std::vector<TString> variables;
  std::vector<TString> unfoldingVariables;
  std::vector<TString> axisTitles;

  std::vector<TString> partonVariables;
  std::vector<TString> partonAxisTitles;

  std::vector<bool> axisInLogScale;

  std::map<TString, TString> extractedSignalFiles;

  std::map<TString, TString> closureScaleFactors;

  std::map<TString, std::map<TString, float>> crossSections;

  TString subDir, baseDir;

  std::vector<float> axisLowValues;
  std::vector<float> axisHighValues;

  const TString lxplusPath = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA";
  const TString localPath = "/Users/georgebakas/Documents/HEP-NTUA_ul";

  std::map<TString, std::map<TString, Double_t>> qcdFitInitValues;

  std::vector<TString> variations;

  bool isUL = true;

  struct correlationMatrix
  {
    Double_t correlations[9];
  };

  std::map<TString, correlationMatrix> correlations;

  std::map<TString, TString> currentlyWorkingDirectory = {
      {"2016_preVFP", "/SR_t02_bM_CR_t02_bM"},
      {"2016_postVFP", "/SR_t02_bM_CR_t02_bM"}, //2016
      {"2017", "/SR_t00_bM_CR_t00_bM"},         //2017
      {"2018", "/SR_t01_bM_CR_t01_bM"}          //2018
  };

  void appendPdfVariations(std::vector<TString> &variations)
  {
    for (int i = 0; i < 100; i++)
    {
      variations.push_back(TString::Format("pdfVariation%i", i));
      correlations.insert(std::pair<TString, correlationMatrix>(TString::Format("pdfVariation%i", i),
                                                                {1., 1., 1.,
                                                                1., 1., 1.,
                                                                1., 1., 1.}));
    }
  }

  void appendScaleVariations(std::vector<TString> &variations)
  {
    for (int i = 1; i < 6; i++)
    {
      variations.push_back(TString::Format("scaleWeight%i", i));
      correlations.insert(std::pair<TString, correlationMatrix>(TString::Format("scaleWeight%i", i),
                                                                {1., 1., 1.,
                                                                1., 1., 1.,
                                                                1., 1., 1.}));
    }
  }

  void appendJESWeights(std::vector<TString> &variations)
  {
    //variations.push_back("boostedSmeared");
    variations.push_back("boostedSmearedUp");
    variations.push_back("boostedShiftedUp");
    variations.push_back("boostedSmearedDown");
    variations.push_back("boostedShiftedDown");
    for (int i = 1; i <= 30; i++)
    {
      variations.push_back(TString::Format("boostedShiftedSrc%iUp", i));
      variations.push_back(TString::Format("boostedShiftedSrc%iDown", i));
    }
  }

  void clearConstants()
  {
    floatConstants.clear();
    luminositiesCR.clear();
    luminositiesSR.clear();
    fitConstants.clear();
    fitConstantsErrors.clear();
    unfoldingVariables.clear();
    variables.clear();
    axisTitles.clear();
    partonVariables.clear();
    partonAxisTitles.clear();

    axisInLogScale.clear();
    extractedSignalFiles.clear();
    closureScaleFactors.clear();
    crossSections.clear();

    variations.clear();

    qcdFitInitValues.clear();
  }

  void initConstants()
  {
    std::cout << "Skata" << std::endl;

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("2.2f");
    clearConstants();
    baseDir = lxplusPath;
    //baseDir = localPath;
    years = {"2016_preVFP", "2016_postVFP", "2017", "2018"};
    /*
    if (!gSystem->AccessPathName(lxplusPath))
    {
      baseDir = lxplusPath;
    }
    else if (!gSystem->AccessPathName(localPath))
    {
      baseDir = localPath;
    }
    else
    {
      std::cout << "Not in lxplus or ipapakri laptop, analysis will not run correctly." << std::endl;
    } */ 
    std::cout<<"initialized everything..."<<std::endl;
    /*
    extractedSignalFiles.insert(std::pair<TString, TString>("2016", baseDir + "/SignalExtraction/results/2016" + currentlyWorkingDirectory["2016"] + "/ExtractedSignal.root"));
    extractedSignalFiles.insert(std::pair<TString, TString>("2017", baseDir + "/SignalExtraction/results/2017/ExtractedSignal.root"));
    extractedSignalFiles.insert(std::pair<TString, TString>("2018", baseDir + "/SignalExtraction/results/2018/ExtractedSignal.root"));

    closureScaleFactors.insert(std::pair<TString, TString>("2016", baseDir + "/ClosureTests/results/2016" + currentlyWorkingDirectory["2016"] + "/reduced/fitted/closureTest_fitResults_2016.root"));
    closureScaleFactors.insert(std::pair<TString, TString>("2017", baseDir + "/ClosureTests/results/2017" + currentlyWorkingDirectory["2017"] + "/reduced/fitted/closureTest_fitResults_2017.root"));
    closureScaleFactors.insert(std::pair<TString, TString>("2018", baseDir + "/ClosureTests/results/2018" + currentlyWorkingDirectory["2018"] + "/reduced/fitted/closureTest_fitResults_2018.root"));
    */
    floatConstants.insert(std::pair<TString, float>("bTagEff2016", 0.62990543));
    floatConstants.insert(std::pair<TString, float>("bTagEff2017", 0.60554176));
    floatConstants.insert(std::pair<TString, float>("bTagEff2018", 0.63387173));

    luminositiesSR.insert(std::pair<TString, float>("2016_preVFP", 19500));
    luminositiesSR.insert(std::pair<TString, float>("2016_postVFP", 16500));
    luminositiesSR.insert(std::pair<TString, float>("2017", 41480));
    luminositiesSR.insert(std::pair<TString, float>("2018", 59830));

    luminositiesCR.insert(std::pair<TString, float>("2016_preVFP", 1134));
    luminositiesCR.insert(std::pair<TString, float>("2016_postVFP", 564));
    luminositiesCR.insert(std::pair<TString, float>("2017", 41480));
    luminositiesCR.insert(std::pair<TString, float>("2018", 59830));
    std::cout<<"initialized everything..."<<std::endl;
    //unfoldingVariables.push_back("leadingJetPt");
    //unfoldingVariables.push_back("subleadingJetPt");
    unfoldingVariables.push_back("jetPt0");
    unfoldingVariables.push_back("jetPt1");
    unfoldingVariables.push_back("yJJ");
    unfoldingVariables.push_back("ptJJ");
    unfoldingVariables.push_back("mJJ");
    //unfoldingVariables.push_back("leadingJetY");
    //unfoldingVariables.push_back("subleadingJetY");
    unfoldingVariables.push_back("jetY0");
    unfoldingVariables.push_back("jetY1");
    //unfoldingVariables.push_back("leadingJetCosThetaStar");
    //unfoldingVariables.push_back("subleadingJetCosThetaStar");
    //unfoldingVariables.push_back("YStar");
    unfoldingVariables.push_back("cosTheta_0");
    unfoldingVariables.push_back("cosTheta_1");
    unfoldingVariables.push_back("chi");
    std::cout<<"initialized everything..."<<std::endl;
    //variables.push_back("leadingJetPt");
    //variables.push_back("subleadingJetPt");
    variables.push_back("jetPt0");
    variables.push_back("jetPt1");
    variables.push_back("yJJ");
    variables.push_back("ptJJ");
    variables.push_back("mJJ");
    //variables.push_back("leadingJetY");
    //variables.push_back("subleadingJetY");
    variables.push_back("jetY0");
    variables.push_back("jetY1");
    //variables.push_back("leadingJetCosThetaStar");
    //variables.push_back("subleadingJetCosThetaStar");
    //variables.push_back("YStar");
    variables.push_back("cosTheta_0");
    variables.push_back("cosTheta_1");
    variables.push_back("chi");
    variables.push_back("leadingJetMassSoftDrop");
    variables.push_back("subleadingJetMassSoftDrop");
    
    partonVariables.push_back("leadingPartonPt");
    partonVariables.push_back("subleadingPartonPt");
    partonVariables.push_back("yJJ");
    partonVariables.push_back("ptTTbarParton");
    partonVariables.push_back("mJJ");
    partonVariables.push_back("leadingPartonY");
    partonVariables.push_back("subleadingPartonY");
    partonVariables.push_back("leadingPartonCosThetaStar");
    partonVariables.push_back("subleadingPartonCosThetaStar");
    partonVariables.push_back("YStar");
    partonVariables.push_back("leadingPartonMass");
    partonVariables.push_back("subleadingPartonMass");

    axisTitles.push_back("Leading jet Pt (GeV)");
    axisTitles.push_back("Subleading jet Pt (GeV)");
    axisTitles.push_back("yJJ");
    axisTitles.push_back("ptJJ (GeV)");
    axisTitles.push_back("mJJ (GeV)");
    axisTitles.push_back("Leading jet Y");
    axisTitles.push_back("Subleading jet Y");
    axisTitles.push_back("leadingPartonCosThetaStar");
    axisTitles.push_back("subleadingPartonCosThetaStar");
    axisTitles.push_back("YStar");
    axisTitles.push_back("Leading jet mass (GeV)");
    axisTitles.push_back("Subleading jet mass (GeV)");

    partonAxisTitles.push_back("Leading parton Pt (GeV)");
    partonAxisTitles.push_back("Subleading parton Pt (GeV)");
    partonAxisTitles.push_back("yJJ");
    partonAxisTitles.push_back("ptJJ (GeV)");
    partonAxisTitles.push_back("mTTbarParton (GeV)");
    partonAxisTitles.push_back("Leading parton Y");
    partonAxisTitles.push_back("Subleading parton Y");
    partonAxisTitles.push_back("leadingPartonCosThetaStar");
    partonAxisTitles.push_back("subleadingPartonCosThetaStar");
    partonAxisTitles.push_back("YStar");
    partonAxisTitles.push_back("Leading parton mass (GeV)");
    partonAxisTitles.push_back("Subleading parton mass (GeV)");

    axisInLogScale.push_back(true);
    axisInLogScale.push_back(true);
    axisInLogScale.push_back(true);
    axisInLogScale.push_back(true);
    axisInLogScale.push_back(true);
    axisInLogScale.push_back(false);
    axisInLogScale.push_back(false);
    axisInLogScale.push_back(true);
    axisInLogScale.push_back(true);

    axisLowValues.push_back(-10);
    axisLowValues.push_back(-10);
    axisLowValues.push_back(-10);
    axisLowValues.push_back(-10);
    axisLowValues.push_back(-10);
    axisLowValues.push_back(0);
    axisLowValues.push_back(0);
    axisLowValues.push_back(-10);
    axisLowValues.push_back(-10);

    axisHighValues.push_back(-10);
    axisHighValues.push_back(-10);
    axisHighValues.push_back(-10);
    axisHighValues.push_back(-10);
    axisHighValues.push_back(-10);
    axisHighValues.push_back(3.5);
    axisHighValues.push_back(3.5);
    axisHighValues.push_back(-10);
    axisHighValues.push_back(-10);

    

    std::map<TString, float> XS2016 = {{"QCD-HT300to500", 315400},
                                      {"QCD-HT500to700", 32260},
                                      {"QCD-HT700to1000", 6830},
                                      {"QCD-HT1000to1500", 1207},
                                      {"QCD-HT1500to2000", 119.1},
                                      {"QCD-HT2000toInf", 25.16},
                                       //{"Mtt-700to1000", 69.64},
                                       //{"Mtt-1000toInf", 16.74},
                                      {"Mtt-700to1000", 80.78},
                                      {"Mtt-1000toInf", 21.43},
                                      {"nominal", 832},
                                      {"Hadronic", 377.96},
                                      {"SemiLeptonic", 365.34},
                                      {"Dilepton", 88.29},
                                      {"DYJetsToQQ", 1208},
                                      {"WJetsToQQ", 3105},
                                      {"ST-tW-top-5f", 38.09},
                                      {"ST-tW-antitop-5f", 38.09},
                                      {"ST-t-channel-top-4f", 35.6},
                                      {"ST-t-channel-antitop-4f", 35.6},
                                      {"ST-t-channel-top-5f", 119.7},
                                      {"ST-t-channel-antitop-5f", 82.52}};

    std::map<TString, float> XS2017 = {{"QCD-HT300to500", 322600},
                                      {"QCD-HT500to700", 29980},
                                      {"QCD-HT700to1000", 6334},
                                      {"QCD-HT1000to1500", 1088},
                                      {"QCD-HT1500to2000", 99.11},
                                      {"QCD-HT2000toInf", 20.23},
                                      {"Mtt-700to1000", 69.64},
                                      {"Mtt-1000toInf", 16.74},
                                      {"nominal", 832},
                                      {"Hadronic", 377.96},
                                      {"SemiLeptonic", 365.34},
                                      {"Dilepton", 88.29},
                                      {"DYJetsToQQ", 1728},
                                      {"WJetsToQQ-HT400to600", 1447},
                                      {"WJetsToQQ-HT600to800", 318.8},
                                      {"ST-tW-top-5f", 34.91},
                                      {"ST-tW-antitop-5f", 34.97},
                                      {"ST-t-channel-top-4f", 115.3},
                                      {"ST-t-channel-antitop-4f", 69.09},
                                      {"ST-t-channel-top-5f", 119.7},
                                      {"ST-t-channel-antitop-5f", 71.74},
                                      {"amcatnlo", 832},
                                      {"madgraph", 832},
                                      {"herwig", 832}};

    std::map<TString, float> XS2018 = {{"QCD-HT300to500", 323400},
                                      {"QCD-HT500to700", 30140},
                                      {"QCD-HT700to1000", 6310},
                                      {"QCD-HT1000to1500", 1094},
                                      {"QCD-HT1500to2000", 99.38},
                                      {"QCD-HT2000toInf", 20.2},
                                      {"Mtt-700to1000", 69.64},
                                      {"Mtt-1000toInf", 16.74},
                                      {"nominal", 832},
                                      {"Hadronic", 377.96},
                                      {"SemiLeptonic", 365.34},
                                      {"Dilepton", 88.29},
                                      {"DYJetsToQQ", 1728},
                                      {"WJetsToQQ-HT400to600", 1447},
                                      {"WJetsToQQ-HT600to800", 318.8},
                                      {"ST-tW-top-5f", 34.91},
                                      {"ST-tW-antitop-5f", 34.97},
                                      {"channel-top-4f", 115.3},
                                      {"channel-antitop-4f", 69.09},
                                      {"channel-antitop-5f", 71.74},
                                      {"channel-top-5f", 119.7},
                                      {"amcatnlo", 832},
                                      {"madgraph", 832},
                                      {"herwig", 832}};

    crossSections.insert(std::pair<TString, std::map<TString, float>>("2016_preVFP", XS2016));
    crossSections.insert(std::pair<TString, std::map<TString, float>>("2016_postVFP", XS2016));
    crossSections.insert(std::pair<TString, std::map<TString, float>>("2017", XS2017));
    crossSections.insert(std::pair<TString, std::map<TString, float>>("2018", XS2018));

    variations = {
        /*"amcatnlo",
        //"madgraph",
        "herwig",
        "hdampdown",
        "hdampup",
        "mtop166",
        "mtop169",
        "mtop171",
        "mtop173",
        "mtop175",
        "mtop178",
        "cp5down",
        "cp5up",
        "isrDefLo",
        "isrDefHi",
        "fsrDefLo",
        "fsrDefHi" ,
        "bTagUp",
        "bTagDown"*/
    };

    //appendJESWeights(variations);
    //appendScaleVariations(variations);
    //appendPdfVariations(variations);

    std::map<TString, Double_t> qcdFitInitValues2016 = {{"b0", 3.0981e-01},
                                                        {"b1", 9.0739e-01},
                                                        {"b2", 2.1589e-07},
                                                        {"b3", 2.2265e-03},
                                                        {"b4", 4.2621e-03}};

    std::map<TString, Double_t> qcdFitInitValues2017 = {{"b0", 3.0981e-01},
                                                        {"b1", 1.9814e+00},
                                                        {"b2", 8.7545e-02},
                                                        {"lb3", 2.2265e-03},
                                                        {"b4", 3.0788e-02}};

    std::map<TString, Double_t> qcdFitInitValues2018 = {{"b0", 2.7840e-01},
                                                        {"b1", 5.6995e-01},
                                                        {"b2", 4.6676e-02},
                                                        {"b3", 1.5639e-03},
                                                        {"b4", 5.8243e-03}};

    qcdFitInitValues.insert(std::pair<TString, std::map<TString, Double_t>>("2016_preVFP", qcdFitInitValues2016));
    qcdFitInitValues.insert(std::pair<TString, std::map<TString, Double_t>>("2016_postVFP", qcdFitInitValues2016));
    qcdFitInitValues.insert(std::pair<TString, std::map<TString, Double_t>>("2017", qcdFitInitValues2017));
    qcdFitInitValues.insert(std::pair<TString, std::map<TString, Double_t>>("2018", qcdFitInitValues2018));
    
    /*
    correlations.insert(std::pair<TString, correlationMatrix>("hdampdown", {1., .0, .0,
                                                                            .0, 1., .0,
                                                                            .0, .0, 1.}));
    correlations.insert(std::pair<TString, correlationMatrix>("hdampup", {1., .0, .0,
                                                                          .0, 1., .0,
                                                                          .0, .0, 1.}));
    correlations.insert(std::pair<TString, correlationMatrix>("cp5down", {1., .0, .0,
                                                                          .0, 1., .0,
                                                                          .0, .0, 1.}));
    correlations.insert(std::pair<TString, correlationMatrix>("cp5up", {1., .0, .0,
                                                                        .0, 1., .0,
                                                                        .0, .0, 1.}));
    correlations.insert(std::pair<TString, correlationMatrix>("isrDefLo", {1., .5, .5,
                                                                           .5, 1., 1.,
                                                                           .5, 1., 1.}));
    correlations.insert(std::pair<TString, correlationMatrix>("isrDefHi", {1., .5, .5,
                                                                           .5, 1., 1.,
                                                                           .5, 1., 1.}));
    correlations.insert(std::pair<TString, correlationMatrix>("fsrDefLo", {1., .5, .5,
                                                                           .5, 1., 1.,
                                                                           .5, 1., 1.}));
    correlations.insert(std::pair<TString, correlationMatrix>("fsrDefHi", {1., .5, .5,
                                                                           .5, 1., 1.,
                                                                           .5, 1., 1.}));
    /*correlations.insert(std::pair<TString, correlationMatrix>("bTagUp", {1., .5, .5,
                                                                         .5, 1., .5,
                                                                         .5, .5, .1}));
    correlations.insert(std::pair<TString, correlationMatrix>("bTagDown", {1., .5, .5,
                                                                           .5, 1., .5,
                                                                           .5, .5, 1.}));
                                                                           */
  }
} // namespace AnalysisConstants

#endif
