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


void initFilesMapping(bool free_eb = true)
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
	
	floatConstants.insert(pair<TString, float>("bTagEff2016", 0.629909));
	floatConstants.insert(pair<TString, float>("bTagEff2017", 0.605622));
	floatConstants.insert(pair<TString, float>("bTagEff2018", 0.633934));


	//these are fit results taken from the simultaneous fit when btagging efficiency eb is let free in the fit
	if(free_eb)
	{
		Nbkg2Constants.insert(pair<TString, float>("Nbkg2016", 3.0802e+03));
		Nbkg2Constants.insert(pair<TString, float>("Nbkg2017", 2.4652e+03));
		Nbkg2Constants.insert(pair<TString, float>("Nbkg2018", 4.4306e+03));

		Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 1.45e+02));
		Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2017_error", 1.32e+02));
		Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 1.77e+02));

		Nbkg0Constants.insert(pair<TString, float>("Nbkg2016", 8.8132e+04));
		Nbkg0Constants.insert(pair<TString, float>("Nbkg2017", 1.5483e+05));
		Nbkg0Constants.insert(pair<TString, float>("Nbkg2018", 1.6842e+05));

		Nbkg0ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 4.48e+02));
		Nbkg0ConstantsErrors.insert(pair<TString, float>("Nbkg2017_error", 5.24e+02));
		Nbkg0ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 5.64e+02));

		//from ABCD method
		NQCD2_reduced.insert(pair<TString, float>("2016",846.365));
		NQCD2_reduced.insert(pair<TString, float>("2017",677.041));
		NQCD2_reduced.insert(pair<TString, float>("2017",455.439));
		NQCD2_reduced.insert(pair<TString, float>("2018",1201.35));

	}
	else
	{
		Nbkg2Constants.insert(pair<TString, float>("Nbkg2016", 3.0552e+03));
		Nbkg2Constants.insert(pair<TString, float>("Nbkg2017", 2.3906e+03));
		Nbkg2Constants.insert(pair<TString, float>("Nbkg2018", 4.3472e+03));

		Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 1.42e+02));
		Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2017_error", 1.45e+02));
		Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 1.78e+02));

		Nbkg0Constants.insert(pair<TString, float>("Nbkg2016", 8.8492e+04));
		Nbkg0Constants.insert(pair<TString, float>("Nbkg2017", 1.5612e+05));
		Nbkg0Constants.insert(pair<TString, float>("Nbkg2018", 1.7002e+05));

		Nbkg0ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 3.70e+02));
		Nbkg0ConstantsErrors.insert(pair<TString, float>("Nbkg2017_error", 3.78e+02));
		Nbkg0ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 4.37e+02));

		//from ABCD method
		NQCD2_reduced.insert(pair<TString, float>("2016",836.08));
		NQCD2_reduced.insert(pair<TString, float>("2017",651.127));
		NQCD2_reduced.insert(pair<TString, float>("2018",1167.64));
	}

	
	
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