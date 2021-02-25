#include "TString.h"

map<TString, map<TString, TString>> files;
map<TString, map<TString, TString>> filesReduced;

map<TString, float> Nbkg2Constants;
map<TString, float> Nbkg2ConstantsErrors;
map<TString, int> variableConstant;
TColor color;

map<TString, float> ttbarSigStrength;
map<TString, float> ttbarSigStrength_noBTagSF;
map<TString, float> ttbarSigStrengthError;
map<TString, float> ttbarSigStrengthError_noBTagSF;

map<TString, float> luminosity;
map<TString, float> luminosityCR;
map<TString, float> deepCSVFloatMap;
map<TString, float> topTaggerCuts;
map<TString, int> triggerSRConst;
map<TString, int> triggerCRConst;


void initFilesMapping()
{
	//thsese will be used either for naming the files for FillHistograms, or either for getting the histo names
	map<TString, TString> files2016_pre = {{"data",  "2016_preVFP/Histo_Data_2016_preVFP_100.root"},
	                                       {"mcSig", "2016_preVFP/Histo_TT_Mtt-700toInf_100.root"},
	                                       {"mcSub", "2016_preVFP/Histo_SubdominantBkgs_100.root"},
	                               	       {"qcd"  , "2016_preVFP/Histo_QCD_HT300toInf_100.root"}};
	
	map<TString, TString> files2016_post = {{"data",  "2016_postVFP/Histo_Data_2016_preVFP_100.root"},
	                                   		{"mcSig", "2016_postVFP/Histo_TT_Mtt-700toInf_100.root"},
	                                   		{"mcSub", "2016_postVFP/Histo_SubdominantBkgs_100.root"},
	                               	  		{"qcd"  , "2016_postVFP/Histo_QCD_HT300toInf_100.root"}};

	map<TString, TString> files2017 = {{"data",  "2017/Histo_Data_2017_100.root"},
	                                   {"mcSig", "2017/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2017/Histo_SubdominantBkgs_100.root"},
	                                   {"qcd"  , "2017/Histo_QCD_HT300toInf_100.root"}};

	map<TString, TString> files2018 = {{"data",  "2018/Histo_Data_2018_100.root"},
	                                   {"mcSig", "2018/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2018/Histo_SubdominantBkgs_100.root"},
	                                   {"qcd"  , "2018/Histo_QCD_HT300toInf_100.root"}};


	files.insert(pair<TString, map<TString, TString>>("2016_preVFP", files2016_pre));
	files.insert(pair<TString, map<TString, TString>>("2016_postVFP", files2016_post));
	files.insert(pair<TString, map<TString, TString>>("2017", files2017));
	files.insert(pair<TString, map<TString, TString>>("2018", files2018));


	map<TString, TString> files2016Reduced_pre = {{"data", "2016_preVFP/Histo_Data_2016_preVFP_100_reduced.root"},
	                                   		  {"mcSig", "2016_preVFP/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016_preVFP/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2016_preVFP/Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2016Reduced_post = {{"data", "2016_postVFP/Histo_Data_2016_postVFP_100_reduced.root"},
	                                   		  {"mcSig", "2016_postVFP/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016_postVFP/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2016_postVFP/Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2017Reduced = {{"data",  "2017/Histo_Data_2017_100_reduced.root"},
	                                   		  {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "20167Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2018Reduced = {{"data",  "2018/Histo_Data_2018_100_reduced.root"},
	                                   	      {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2018/Histo_QCD_HT300toInf_100_reduced.root"}};

	filesReduced.insert(pair<TString, map<TString, TString>>("2016_preVFP", files2016Reduced_pre));
	filesReduced.insert(pair<TString, map<TString, TString>>("2016_postVFP", files2016Reduced_post));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017", files2017Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018", files2018Reduced));

	luminosity.insert(pair<TString, float>("2016_preVFP", 19500));
	luminosity.insert(pair<TString, float>("2016_postVFP", 16500));
	luminosity.insert(pair<TString, float>("2017",41530));
	luminosity.insert(pair<TString, float>("2018",59740));

	luminosityCR.insert(pair<TString, float>("2016_preVFP", 1134));
	luminosityCR.insert(pair<TString, float>("2016_postVFP", 564));
	luminosityCR.insert(pair<TString, float>("2017",41530));
	luminosityCR.insert(pair<TString, float>("2018",59740));

	topTaggerCuts.insert(pair<TString, float>("2016_preVFP", 0.2));
	topTaggerCuts.insert(pair<TString, float>("2016_postVFP", 0.2));
	topTaggerCuts.insert(pair<TString, float>("2017", 0.0));
	topTaggerCuts.insert(pair<TString, float>("2018", 0.1));

	deepCSVFloatMap.insert(pair<TString, float>("2016_preVFP", 0.6321));
	deepCSVFloatMap.insert(pair<TString, float>("2016_postVFP", 0.6321));
	deepCSVFloatMap.insert(pair<TString, float>("2017",0.4506));
	deepCSVFloatMap.insert(pair<TString, float>("2018",0.4168));

	triggerSRConst.insert(pair<TString, int>("2016_preVFP",2));
	triggerSRConst.insert(pair<TString, int>("2016_postVFP",2));
	triggerSRConst.insert(pair<TString, int>("2017",5));
	triggerSRConst.insert(pair<TString, int>("2018",5));

	triggerCRConst.insert(pair<TString, int>("2016_preVFP",4));
	triggerCRConst.insert(pair<TString, int>("2016_postVFP",4));
	triggerCRConst.insert(pair<TString, int>("2017",5));
	triggerCRConst.insert(pair<TString, int>("2018",5));

	variableConstant.insert(pair<TString, int>("mJJ",  0));
	variableConstant.insert(pair<TString, int>("ptJJ", 1));
	variableConstant.insert(pair<TString, int>("yJJ",  2));
	variableConstant.insert(pair<TString, int>("jetPt0", 3));
	variableConstant.insert(pair<TString, int>("jetPt1", 4));
	variableConstant.insert(pair<TString, int>("jetY0", 5));
	variableConstant.insert(pair<TString, int>("jetY1", 6));

}
