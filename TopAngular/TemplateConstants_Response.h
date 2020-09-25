#include "TString.h"

map<TString, map<TString, TString>> files;

map<TString, float> floatConstants;
map<TString, float> topTaggerConstants;
map<TString, int> variableConstant;
map<TString, int> variableConstantParton;
map<TString, int> variableConstantParticle;
map<TString, float> luminosity;
map<TString, TString> eospath;
TColor color;
//std::vector< std::vector <Float_t> > BND;

void initFilesMapping()
{

	map<TString, TString> files2016 = {{"700-1000", "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"},
	                                   {"1000-Inf", "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"},
	                               	   {"TTNominal", "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"},
																		 {"TTHadronic_0", "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
																		 {"TTSemiLeptonic_0", "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
																		 {"TTTo2L2Nu_0", "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root"}};

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

	eospath.insert(pair<TString,TString>("2016","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/"));
	eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2017/Signal/"));
	eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Signal/"));

	files.insert(pair<TString, map<TString, TString>>("2016", files2016));
	files.insert(pair<TString, map<TString, TString>>("2017", files2017));
	files.insert(pair<TString, map<TString, TString>>("2018", files2018));


	//these are btagging Working points for each year Medium WP
	floatConstants.insert(pair<TString, float>("btagWP2016", 0.6321));
	floatConstants.insert(pair<TString, float>("btagWP2017", 0.4941));
	floatConstants.insert(pair<TString, float>("btagWP2018", 0.4184));

	//these are btagging Working points for each year Loose WP
	//floatConstants.insert(pair<TString, float>("btagWP2016", 0.2217));
	//floatConstants.insert(pair<TString, float>("btagWP2017", 0.1522));
	//floatConstants.insert(pair<TString, float>("btagWP2018", 0.1241));

	topTaggerConstants.insert(pair<TString, float>("topTagger2016", 0.2));
	topTaggerConstants.insert(pair<TString, float>("topTagger2017", 0.0));
	topTaggerConstants.insert(pair<TString, float>("topTagger2018", 0.1));

	variableConstant.insert(pair<TString, int>("mJJ",  0));
	variableConstant.insert(pair<TString, int>("ptJJ", 1));
	variableConstant.insert(pair<TString, int>("yJJ",  2));
	variableConstant.insert(pair<TString, int>("jetPt0", 3));
	variableConstant.insert(pair<TString, int>("jetPt1", 4));
	variableConstant.insert(pair<TString, int>("jetY0", 5));
	variableConstant.insert(pair<TString, int>("jetY1", 6));

	variableConstantParton.insert(pair<TString, int>("mTTbarParton",  0));
	variableConstantParton.insert(pair<TString, int>("ptTTbarParton", 1));
	variableConstantParton.insert(pair<TString, int>("yTTbarParton",  2));
	variableConstantParton.insert(pair<TString, int>("ptTopParton0", 3));
	variableConstantParton.insert(pair<TString, int>("ptTopParton1", 4));
	variableConstantParton.insert(pair<TString, int>("yTopParton0", 5));
	variableConstantParton.insert(pair<TString, int>("yTopParton1", 6));

	variableConstantParticle.insert(pair<TString, int>("mJJGen",  0));
	variableConstantParticle.insert(pair<TString, int>("ptJJGen", 1));
	variableConstantParticle.insert(pair<TString, int>("yJJGen",  2));
	variableConstantParticle.insert(pair<TString, int>("genjetPt0", 3));
	variableConstantParticle.insert(pair<TString, int>("genjetPt1", 4));
	variableConstantParticle.insert(pair<TString, int>("genjetY0", 5));
	variableConstantParticle.insert(pair<TString, int>("genjetY1", 6));

	luminosity.insert(pair<TString,float>("luminosity2016", 35920));
	luminosity.insert(pair<TString,float>("luminosity2017", 41530));
	luminosity.insert(pair<TString,float>("luminosity2018", 59740));
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
