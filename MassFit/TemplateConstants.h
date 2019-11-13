#include "TString.h"

map<TString, map<TString, TString>> files;

map<TString, float> floatConstants;

TColor color;

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
	
	floatConstants.insert(pair<TString, float>("bTagEff2016", 0.629909));
	floatConstants.insert(pair<TString, float>("bTagEff2017", 0.605622));
	floatConstants.insert(pair<TString, float>("bTagEff2018", 0.633934));
}

void initNqcdVsEb(TString year)
{

	const int n = 7;
	if(year.EqualTo("2016"))
	{
		float eb[n]= {0.4, 0.5, 0.56029, 0.6, 0.629909, 0.7, 0.8};
		
		float error_x[n] = {0,0, 1.17e-02, 0,0, 0, 0};
		float NQCD[3][n] = {
				{8.1973e+04, 8.5459e+04, 8.7019e+04, 8.7826e+04,8.832e+04 ,8.9197e+04, 9.9892e+04}, //0btag
				{2.6618e+04, 2.8082e+04, 2.8973e+04, 2.9865e+04,3.0542e+04 ,3.2121e+04, 3.4154e+04}, //1btag
				{4.1645e+03,3.2985e+03, 2.9980e+03,2.8784e+03,2.84e+03 ,2.8091e+03, 2.9322e+03} //2btag
				};

		float error_y[3][n] = {
				{4.08e+02, 3.36e+02, 4.15e+02, 3.15e+02, 3.13e+02,3.10e+02, 3.08e+02}, //0btag
				{3.80e+02, 3.20e+02, 3.93e+02, 2.76e+02,2.62e+02 ,2.36e+02, 2.13e+02}, //1btag
				{1.36e+02, 1.37e+02, 1.43e+02, 1.38e+02,1.55e+02 ,1.41e+02, 1.45e+02} //2btag
				};
	}
	else if(year.EqualTo("2017"))
	{
		float eb[n] = {0.4, 0.482, 0.5, 0.6, 0.605622, 0.7, 0.8};
		float error_x[n] = {0,1.08e-02,0,0,0,0,0};

		float NQCD[3][n] = {
				{1.4919e+05, 1.5273e+05, 1.5334e+05, 1.5579e+05, 1.5591e+05, 1.5699e+05, 1.5759e+05}, //0btag
				{3.3164e+04, 3.4635e+04, 3.5087e+04, 3.7612e+04, 3.7752e+04, 3.9268e+04, 4.0439e+04}, //1btag
				{2.8696e+03, 2.3896e+03, 2.3186e+03, 2.1258e+03, 2.1182e+03 , 1.5872e+03, 462.59} //2btag
				};

		float error_y[3][n]= {
				{475, 571, 424, 409, 408, 404, 404}, //0btag
				{387, 431, 347, 284, 283, 217, 216}, //1btag
				{1.27e+02, 1.36e+02, 1.28e+02, 1.33e+02, 1.33e+02 , 0.730e+02, 0.07e+02} //2btag
				};
	}
	else
	{
		float eb[n] = {0.4, 0.5, 0.51255, 0.6, 0.633934, 0.7, 0.8};
		float error_x[n] = {0, 0,9.8e-03, 0,0, 0, 0};

		float NQCD[3][n] = {
				{1.5893e+05, 1.6500e+05, 1.6558e+05, 1.6885e+05, 1.6973e+05, 1.7098e+05, 1.7192e+05}, //0btag
				{4.2713e+04, 4.5549e+04, 4.5743e+04, 4.8995e+04, 5.0264e+04, 5.2599e+04, 5.4821e+04}, //1btag
				{5.4029e+03, 4.2967e+03, 4.2092e+03, 3.8442e+03, 3.7944e+03 , 3.8195e+03, 2.7492e+03} //2btag
				};

		float error_y[3][n]= {
				{544,456,668,429,426,419,421}, //0btag
				{498,334,278,343,325,383,250}, //1btag
				{173,170,142,178,181,120,101}//2btag
				};

	}

}