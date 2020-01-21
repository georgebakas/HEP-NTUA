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
	

	map<TString, TString> filesTT2016 = {{"700-1000", "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"},
	                                   {"1000-Inf", "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"}};
																		 
	map<TString, TString> filesTT2017 = {{"700-1000", "TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8.root"},
	                                   {"1000-Inf", "TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root"}};
	
	map<TString, TString> filesTT2018 = {{"700-1000", "TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
	                                   {"1000-Inf", "TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root"}};
																		 
	filesTT.insert(pair<TString, map<TString, TString>>("2016", filesTT2016));
	filesTT.insert(pair<TString, map<TString, TString>>("2017", filesTT2017));
	filesTT.insert(pair<TString, map<TString, TString>>("2018", filesTT2018));


	map<TString, TString> filesQCD2016 = {{"300-500", "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                                   	  {"500-700", "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                                   	  {"700-1000", "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                                      {"1000-1500", "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                                      {"1500-2000", "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                                  	  {"2000-Inf", "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"}};
																		 
	map<TString, TString> filesQCD2017 = {{"300-500", "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root"},
	                                   	  {"500-700", "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.root"},
	                                   	  {"700-1000", "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.root"},
	                                      {"1000-1500", "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.root"},
	                                      {"1500-2000", "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root"},
	                                  	  {"2000-Inf", "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.root"}};
	
	map<TString, TString> filesQCD2018 = {{"300-500", "QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                                   	  {"500-700", "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                                   	  {"700-1000", "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                                      {"1000-1500", "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                                      {"1500-2000", "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                                  	  {"2000-Inf", "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root"}};
																		 
	filesQCD.insert(pair<TString, map<TString, TString>>("2016", filesQCD2016));
	filesQCD.insert(pair<TString, map<TString, TString>>("2017", filesQCD2017));
	filesQCD.insert(pair<TString, map<TString, TString>>("2018", filesQCD2018));


	filesSignal.insert(pair<TString,TString>("2016","SignalOutput_AllRegions_0.20_deepCSV_2016.root"));
	filesSignal.insert(pair<TString,TString>("2017","SignalOutput_AllRegions_0.00_deepCSV_2017.root"));
	filesSignal.insert(pair<TString,TString>("2018","SignalOutput_AllRegions_0.10_deepCSV_2018.root"));
	filesSignal.insert(pair<TString,TString>("2016_Loose","SignalOutput_AllRegions_0.20_deepCSVLoose_2016.root"));
	filesSignal.insert(pair<TString,TString>("2017_Loose","SignalOutput_AllRegions_0.00_deepCSVLoose_2017.root"));
	filesSignal.insert(pair<TString,TString>("2018_Loose","SignalOutput_AllRegions_0.10_deepCSVLoose_2018.root"));

	filesBkg.insert(pair<TString,TString>("2016","BkgOutput_AllRegions_0.20_deepCSV_2016.root"));
	filesBkg.insert(pair<TString,TString>("2017","BkgOutput_AllRegions_0.00_deepCSV_2017.root"));
	filesBkg.insert(pair<TString,TString>("2018","BkgOutput_AllRegions_0.10_deepCSV_2018.root"));
	filesBkg.insert(pair<TString,TString>("2016_Loose","BkgOutput_AllRegions_0.20_deepCSVLoose_2016.root"));
	filesBkg.insert(pair<TString,TString>("2017_Loose","BkgOutput_AllRegions_0.00_deepCSVLoose_2017.root"));
	filesBkg.insert(pair<TString,TString>("2018_Loose","BkgOutput_AllRegions_0.10_deepCSVLoose_2018.root"));

	triggerConstant.insert(pair<TString,int>("2016", 2));
	triggerConstant.insert(pair<TString,int>("2017", 5));
	triggerConstant.insert(pair<TString,int>("2018", 5));

	luminosity.insert(pair<TString, float>("2016", 35920));
	luminosity.insert(pair<TString, float>("2017", 41530));
	luminosity.insert(pair<TString, float>("2018", 59740));

	//these are Medium WP's
	deepCSVWP.insert(pair<TString, float>("2016", 0.6321));
	deepCSVWP.insert(pair<TString, float>("2017", 0.4941));
	deepCSVWP.insert(pair<TString, float>("2018", 0.4184));

	//these are Loose WP's
	//deepCSVWP.insert(pair<TString, float>("2016", 0.2217));
	//deepCSVWP.insert(pair<TString, float>("2017", 0.1522));
	//deepCSVWP.insert(pair<TString, float>("2018", 0.1241));

	tTaggerSel.insert(pair<TString, float>("2016", 0.2));
	tTaggerSel.insert(pair<TString, float>("2017", 0.0));
	tTaggerSel.insert(pair<TString, float>("2018", 0.1));


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
