/*
	For Unfolding we need the files:
	1. Signal Extracted S_i(xReco) from SignalExtraction.cpp in MassFit folder 
		MassFit/year/FiducialMeasurements/ABCDMethod/free_eb/SignalHistograms_ABCDMethod_freeEb.root
		we need to rebin the histogram taken from the S_i because it has bins same as parton
	2. Acceptance from files in ../ResponseMatrices/year/UnequalBins/ResponsesEfficiency_year.root more bins for acc
	3. Response matrix from ../ResponseMatrices/year/UnequalBins/ResponsesEfficiency_year.root 
	4. Efficiency from same files but with fewer bins
	5. Lumi taken from this file for every year

*/

#include <TString.h>
map<TString, map<TString, TString>> files;
map<TString, map<TString, TString>> filesReduced;

map<TString, float> floatConstants;

map<TString, float> Nbkg2Constants;
map<TString, float> Nbkg0Constants;
map<TString, float> Nbkg2ConstantsErrors;
map<TString, float> Nbkg0ConstantsErrors;
map<TString, int> variableConstant;  
map<TString, float> NQCD2_reduced;
TColor color;

map<TString, float> luminosity;
map<TString, TString> eospath;
//std::vector< std::vector <Float_t> > BND;




void initFilesMapping(bool free_eb = true)
{
	
	floatConstants.insert(pair<TString, float>("bTagEff2016", 0.629909));
	floatConstants.insert(pair<TString, float>("bTagEff2017", 0.605622));
	floatConstants.insert(pair<TString, float>("bTagEff2018", 0.633934));


	map<TString, TString> files2016 = {{"700-1000", "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"},
	                                {"1000-Inf", "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"},
	                               	{"TTNominal", "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"}};

	map<TString, TString> files2017 = {{"700-1000", "TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8.root"},
	                                {"1000-Inf", "TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root"},
	                                {"TTHadronic_0", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	{"TTSemiLeptonic_0", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	{"TTTo2L2Nu_0", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};
	
	map<TString, TString> files2018 = {{"700-1000", "TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
	                                {"1000-Inf", "TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root"},
	                                {"TTHadronic_0", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	{"TTSemiLeptonic_0", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	{"TTTo2L2Nu_0", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};
										
	eospath.insert(pair<TString,TString>("2016","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_preVFP/Signal/"));
	eospath.insert(pair<TString,TString>("2016","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_postVFP/Signal/"));
	eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Signal/"));
	eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2018/Signal/"));   

	files.insert(pair<TString, map<TString, TString>>("2016_preVFP", files2016));
	files.insert(pair<TString, map<TString, TString>>("2016_postVFP", files2016));
	files.insert(pair<TString, map<TString, TString>>("2017", files2017));
	files.insert(pair<TString, map<TString, TString>>("2018", files2018));

	
	luminosity.insert(pair<TString, float>("2016_preVFP",19500));
	luminosity.insert(pair<TString, float>("2016_postVFP",16500));
	luminosity.insert(pair<TString, float>("2017",41480));
	luminosity.insert(pair<TString, float>("2018",59830));

	variableConstant.insert(pair<TString, int>("mJJ",  0));
	variableConstant.insert(pair<TString, int>("ptJJ", 1));
	variableConstant.insert(pair<TString, int>("yJJ",  2));
	variableConstant.insert(pair<TString, int>("jetPt0", 3));
	variableConstant.insert(pair<TString, int>("jetPt1", 4));
	variableConstant.insert(pair<TString, int>("jetY0", 5));
	variableConstant.insert(pair<TString, int>("jetY1", 6));
	/*
	BND[variableConstant["mJJ"]].push_back({1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}); //mjj
	BND[variableConstant["ptJJ"]].push_back({0,60,150,300,450,600,750,950,1100,1300}); //ptjj
	BND[variableConstant["yJJ"]].push_back({-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}); //yjj
	BND[variableConstant["jetPt0"]].push_back({400,450,500,570,650,750,850,950,1100,1300,1500}); //jetPt0	
	BND[variableConstant["jetPt1"]].push_back({400,450,500,570,650,750,850,950,1100,1300,1500}); //jetPt1
	BND[variableConstant["jetY0"]].push_back({0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}); //jetY0
	BND[variableConstant["jetY1"]].push_back({0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}); //jetY1
	*/
}
