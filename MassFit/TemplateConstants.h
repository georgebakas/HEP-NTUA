#include "TString.h"

map<TString, map<TString, TString>> files;
map<TString, map<TString, TString>> filesReduced;

map<TString, float> Nbkg2Constants;
map<TString, float> Nbkg2ConstantsErrors;
map<TString, int> variableConstant;
TColor color;

map<TString, float> luminosity;
map<TString, float> floatConstants;
map<TString, float> ttbarSigStrength;
map<TString, float> ttbarSigStrength_noBTagSF;
map<TString, float> ttbarSigStrengthError;
map<TString, float> ttbarSigStrengthError_noBTagSF;


void initFilesMapping()
{
	map<TString, TString> files2016 = {{"data",  "2016/Histo_Data_2016_100.root"},
	                                   {"mcSig", "2016/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2016/Histo_SubdominantBkgs_100.root"}};

	map<TString, TString> files2017 = {{"data",  "2017/Histo_Data_2017_100.root"},
	                                   {"mcSig", "2017/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2017/Histo_SubdominantBkgs_100.root"}};

	map<TString, TString> files2018 = {{"data",  "2018/Histo_Data_2018_100.root"},
	                                   {"mcSig", "2018/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2018/Histo_SubdominantBkgs_100.root"}};

	files.insert(pair<TString, map<TString, TString>>("2016", files2016));
	files.insert(pair<TString, map<TString, TString>>("2017", files2017));
	files.insert(pair<TString, map<TString, TString>>("2018", files2018));

	map<TString, TString> files2016Reduced = {{"data",  "2016/Histo_Data_2016_100_reduced.root"},
	                                   		  {"mcSig", "2016/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016/Histo_SubdominantBkgs_100_reduced.root"}};

	map<TString, TString> files2017Reduced = {{"data",  "2017/Histo_Data_2017_100_reduced.root"},
	                                   		  {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced.root"}};

	map<TString, TString> files2018Reduced = {{"data",  "2018/Histo_Data_2018_100_reduced.root"},
	                                   	      {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced.root"}};

	filesReduced.insert(pair<TString, map<TString, TString>>("2016", files2016Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017", files2017Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018", files2018Reduced));


	luminosity.insert(pair<TString, float>("2016",35920));
	luminosity.insert(pair<TString, float>("2017",41530));
	luminosity.insert(pair<TString, float>("2018",59740));

	variableConstant.insert(pair<TString, int>("mJJ",  0));
	variableConstant.insert(pair<TString, int>("ptJJ", 1));
	variableConstant.insert(pair<TString, int>("yJJ",  2));
	variableConstant.insert(pair<TString, int>("jetPt0", 3));
	variableConstant.insert(pair<TString, int>("jetPt1", 4));
	variableConstant.insert(pair<TString, int>("jetY0", 5));
	variableConstant.insert(pair<TString, int>("jetY1", 6));

	//these are fit results from mixed medium wp and loose wp
	//Nbkg2Constants.insert(pair<TString, float>("Nbkg2016", 2.9886e+03));
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2016", 2.9890e+03));
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2017", 2.1662e+03));
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2018", 4.9706e+03));
	//Nbkg2Constants.insert(pair<TString, float>("Nbkg2018", 4.7747e+03));


	//Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 1.73e+02));
	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 1.74e+02));
	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2017_error", 3.11e+02));
	//Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 3.04e+02));
	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 2.95e+02));

	floatConstants.insert(pair<TString, float>("bTagEff2016", 0.629909));
	floatConstants.insert(pair<TString, float>("bTagEff2017", 0.605622));
	floatConstants.insert(pair<TString, float>("bTagEff2018", 0.633934));

	ttbarSigStrength.insert(pair<TString, float>("2016", 0.686668));
	ttbarSigStrength.insert(pair<TString, float>("2017", 0.644361));
	ttbarSigStrength.insert(pair<TString, float>("2018", 0.686214));

	ttbarSigStrength_noBTagSF.insert(pair<TString, float>("2016", 0.671244));
	ttbarSigStrength_noBTagSF.insert(pair<TString, float>("2017", 0.553099));
	ttbarSigStrength_noBTagSF.insert(pair<TString, float>("2018", 0.615816));

	ttbarSigStrengthError.insert(pair<TString, float>("2016", 0.0263103));
	ttbarSigStrengthError.insert(pair<TString, float>("2017", 0.023851));
	ttbarSigStrengthError.insert(pair<TString, float>("2018", 0.019771));

	ttbarSigStrengthError_noBTagSF.insert(pair<TString, float>("2016", 0.0252439));
	ttbarSigStrengthError_noBTagSF.insert(pair<TString, float>("2017", 0.0198563));
	ttbarSigStrengthError_noBTagSF.insert(pair<TString, float>("2018", 0.017298));


}
