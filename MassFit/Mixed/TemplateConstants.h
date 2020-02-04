#include "TString.h"

map<TString, map<TString, TString>> files;
map<TString, map<TString, TString>> filesReduced;

map<TString, float> Nbkg2Constants;
map<TString, float> Nbkg2ConstantsErrors;
map<TString, int> variableConstant;  
TColor color;

map<TString, float> luminosity;
map<TString, float> floatConstants	;


void initFilesMapping()
{
	map<TString, TString> files2016 = {{"data",  "2016/Histo_Data_2016_100.root"},
	                                   {"mcSig", "2016/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2016/Histo_SubdominantBkgs_100.root"}};

	map<TString, TString> files2016Loose = {{"data",  "2016/Histo_Data_2016_100_Loose.root"},
	                                   		{"mcSig", "2016/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   		{"mcSub", "2016/Histo_SubdominantBkgs_100_Loose.root"}};	                                   
																		 
	map<TString, TString> files2017 = {{"data",  "2017/Histo_Data_2017_100.root"},
	                                   {"mcSig", "2017/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2017/Histo_SubdominantBkgs_100.root"}};
    
    map<TString, TString> files2017Loose = {{"data",  "2017/Histo_Data_2017_100_Loose.root"},
	                                   		{"mcSig", "2017/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   		{"mcSub", "2017/Histo_SubdominantBkgs_100_Loose.root"}};		                              
	
	map<TString, TString> files2018 = {{"data",  "2018/Histo_Data_2018_100.root"},
	                                   {"mcSig", "2018/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2018/Histo_SubdominantBkgs_100.root"}};
																		 
	map<TString, TString> files2018Loose = {{"data",  "2018/Histo_Data_2018_100_Loose.root"},
	                                   		{"mcSig", "2018/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   		{"mcSub", "2018/Histo_SubdominantBkgs_100_Loose.root"}};	
	                                   																			 
	files.insert(pair<TString, map<TString, TString>>("2016", files2016));
	files.insert(pair<TString, map<TString, TString>>("2017", files2017));
	files.insert(pair<TString, map<TString, TString>>("2018", files2018));
	files.insert(pair<TString, map<TString, TString>>("2016_Loose", files2016Loose));
	files.insert(pair<TString, map<TString, TString>>("2017_Loose", files2017Loose));
	files.insert(pair<TString, map<TString, TString>>("2018_Loose", files2018Loose));


	map<TString, TString> files2016Reduced = {{"data",  "2016/Histo_Data_2016_100_reduced.root"},
	                                   		  {"mcSig", "2016/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016/Histo_SubdominantBkgs_100_reduced.root"}};
																		 
	map<TString, TString> files2017Reduced = {{"data",  "2017/Histo_Data_2017_100_reduced.root"},
	                                   		  {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced.root"}};                                  		  
	
	map<TString, TString> files2018Reduced = {{"data",  "2018/Histo_Data_2018_100_reduced.root"},
	                                   	      {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced.root"}};

	map<TString, TString> files2016ReducedLoose = {{"data",  "2016/Histo_Data_2016_100_reduced_Loose.root"},
	                                   		       {"mcSig", "2016/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		       {"mcSub", "2016/Histo_SubdominantBkgs_100_reduced_Loose.root"}};
																		 
	map<TString, TString> files2017ReducedLoose = {{"data",  "2017/Histo_Data_2017_100_reduced_Loose.root"},
	                                   		  	   {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		  	   {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced_Loose.root"}};                                  		  
	
	map<TString, TString> files2018ReducedLoose = {{"data",  "2018/Histo_Data_2018_100_reduced_Loose.root"},
	                                   	           {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		       {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced_Loose.root"}};	                                   		  
																		 
	filesReduced.insert(pair<TString, map<TString, TString>>("2016", files2016Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017", files2017Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018", files2018Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2016_Loose", files2016ReducedLoose));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017_Loose", files2017ReducedLoose));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018_Loose", files2018ReducedLoose));
	
	
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
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2016", 2.7202e+03));
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2017", 2.6557e+03));
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2018", 4.7447e+03));

	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 1.91e+02));
	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2017_error", 2.55e+02));
	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 3.27e+02));

	floatConstants.insert(pair<TString, float>("bTagEff2016", 0.629909));
	floatConstants.insert(pair<TString, float>("bTagEff2017", 0.605622));
	floatConstants.insert(pair<TString, float>("bTagEff2018", 0.633934));
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