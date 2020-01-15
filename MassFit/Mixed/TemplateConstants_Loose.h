#include "TString.h"

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
//std::vector< std::vector <Float_t> > BND;


void initFilesMapping()
{

	map<TString, TString> files2016Loose = {{"data",  "2016_Loose/Histo_Data_2016_Loose_100.root"},
	                                   		{"mcSig", "2016_Loose/Histo_TT_Mtt-700toInf_100.root"},
	                                   		{"mcSub", "2016_Loose/Histo_SubdominantBkgs_100.root"}};	                                   
	/*															 
	map<TString, TString> files2017 = {{"data",  "2017_Loose/Histo_Data_2017_100.root"},
	                                   {"mcSig", "2017_Loose/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2017_Loose/Histo_SubdominantBkgs_100.root"}};
	
	map<TString, TString> files2018 = {{"data",  "2018_Loose/Histo_Data_2018_100.root"},
	                                   {"mcSig", "2018_Loose/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2018_Loose/Histo_SubdominantBkgs_100.root"}};
	*/									 
	files.insert(pair<TString, map<TString, TString>>("2016_Loose", files2016Loose));
	//files.insert(pair<TString, map<TString, TString>>("2017_Loose", files2017));
	//files.insert(pair<TString, map<TString, TString>>("2018_Loose", files2018));


	map<TString, TString> files2016Reduced = {{"data",  "2016_Loose/Histo_Data_2016_Loose_100_reduced.root"},
	                                   		  {"mcSig", "2016_Loose/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016_Loose/Histo_SubdominantBkgs_100_reduced.root"}};
	
	/*																 
	map<TString, TString> files2017Reduced = {{"data",  "2017_Loose/Histo_Data_2017_100_reduced.root"},
	                                   		  {"mcSig", "2017_Loose/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2017_Loose/Histo_SubdominantBkgs_100_reduced.root"}};                                  		  
	
	map<TString, TString> files2018Reduced = {{"data",  "2018_Loose/Histo_Data_2018_100_reduced.root"},
	                                   	      {"mcSig", "2018_Loose/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2018_Loose/Histo_SubdominantBkgs_100_reduced.root"}};
	*/																 
	filesReduced.insert(pair<TString, map<TString, TString>>("2016_Loose", files2016Reduced));
	//filesReduced.insert(pair<TString, map<TString, TString>>("2017_Loose", files2017Reduced));
	//filesReduced.insert(pair<TString, map<TString, TString>>("2018_Loose", files2018Reduced));
	
	//floatConstants.insert(pair<TString, float>("bTagEff2016_Loose", 0.629909));
	//floatConstants.insert(pair<TString, float>("bTagEff2017_Loose", 0.605622));
	//floatConstants.insert(pair<TString, float>("bTagEff2018_Loose", 0.633934));


	//these are fit results taken from the simultaneous fit when btagging efficiency eb is let free in the fit
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2016_Loose", 2.2646e+04));
	//Nbkg2Constants.insert(pair<TString, float>("Nbkg2017_Loose", 2.4407e+03));
	//Nbkg2Constants.insert(pair<TString, float>("Nbkg2018_Loose", 4.4094e+03));

	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2016_Loose_error", 6.85e+02));
	//Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2017_Loose_error", 1.31e+02));
	//Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2018_Loose_error", 1.77e+02));


	//from ABCD method
	//NQCD2_reduced.insert(pair<TString, float>("2016",819.074));
	//NQCD2_reduced.insert(pair<TString, float>("2017",659.932));
	//NQCD2_reduced.insert(pair<TString, float>("2018",1176.09));
	
	luminosity.insert(pair<TString, float>("2016_Loose",35920));
	//luminosity.insert(pair<TString, float>("2017_Loose",41530));
	//luminosity.insert(pair<TString, float>("2018_Loose",59740));

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