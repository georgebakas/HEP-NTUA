#include "TString.h"

map<TString, map<TString, TString>> filesTT;
map<TString, map<TString, TString>> filesQCD;

map<TString, int> variableConstant;  
map<TString, int> triggerConstant; 
map<TString, float> luminosity;
map<TString, float> deepCSVWP;
map<TString, float> tTaggerSel;

map<TString, TString> filesSignal;
map<TString, TString> filesBkg;

void initFilesMapping()
{


	filesSignal.insert(pair<TString,TString>("2016","../MassFit/2016/Histo_TT_NominalMC_100.root"));
	filesSignal.insert(pair<TString,TString>("2017","../MassFit/2017/Histo_TT_NominalMC_100.root"));
	filesSignal.insert(pair<TString,TString>("2018","../MassFit/2018/Histo_TT_NominalMC_100.root"));

	filesBkg.insert(pair<TString,TString>("2016","../MassFit/2016/Histo_QCD_HT300toInf_100.root"));
	filesBkg.insert(pair<TString,TString>("2017","../MassFit/2017/Histo_QCD_HT300toInf_100.root"));
	filesBkg.insert(pair<TString,TString>("2018","../MassFit/2018/Histo_QCD_HT300toInf_100.root"));


	luminosity.insert(pair<TString, float>("2016", 35920));
	luminosity.insert(pair<TString, float>("2017", 41530));
	luminosity.insert(pair<TString, float>("2018", 59740));
	
	//these are Medium WP's for Btagging
	deepCSVWP.insert(pair<TString, float>("2016", 0.6321));
	deepCSVWP.insert(pair<TString, float>("2017", 0.4941));
	deepCSVWP.insert(pair<TString, float>("2018", 0.4184));

	tTaggerSel.insert(pair<TString, float>("2016", 0.2));
	tTaggerSel.insert(pair<TString, float>("2017", 0.0));
	tTaggerSel.insert(pair<TString, float>("2018", 0.1));

}
