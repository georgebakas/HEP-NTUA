#include "TString.h"

map<TString, map<TString, TString>> files;
map<TString, map<TString, float>> XSECAll;

map<TString, float> floatConstants;
map<TString, float> topTaggerConstants;
map<TString, float> luminosity;
map<TString, TString> eospath;

void initFilesMapping(bool isTTbar)
{
	if(!isTTbar){
		eospath.insert(pair<TString,TString>("2016","/eos/cms/store/user/gbakas/ZprimeToTT/mc-2016_btag/"));
		eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ZprimeToTT/mc-2017_btag/"));
		eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ZprimeToTT/mc-2018_btag/"));


	map<TString, float> XSEC2016 = {{"M-1000", 0.5392},
																	{"M-1500", 0.07238},
																	{"M-2000", 0.01375},
																	{"M-2500", 0.003},
																	{"M-3000", 0.0007569},
																	{"M-3500", 0.0001979},
																	{"M-4000", 5.484e-05},
																	{"M-5000", 5.043e-06}};

	map<TString, float> XSEC2017 = {{"M1000", 0.5767},
																	{"M1500", 0.07772},
																	{"M2000", 0.01463},
																	{"M2500", 0.003297},
																	{"M3000", 0.0008249},
																	{"M3500", 0.0002226},
																	{"M4000", 6.209e-05},
																	{"M5000", 5.043e-06}};


	map<TString, float> XSEC2018 = {{"M1000", 0.578},
																	{"M1500", 0.07794},
																	{"M2000", 0.01459},
																	{"M2500", 0.003286},
																	{"M3000", 0.0008265},
																	{"M3500", 0.0002228},
																	{"M4000", 6.214e-05},
																	{"M5000", 5.053e-06}};

	XSECAll.insert(pair<TString, map<TString, float>>("2016", XSEC2016));
	XSECAll.insert(pair<TString, map<TString, float>>("2017", XSEC2017));
	XSECAll.insert(pair<TString, map<TString, float>>("2018", XSEC2018));
	}
	else{
		eospath.insert(pair<TString,TString>("2016","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/"));
		eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2017/Signal/"));
		eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Signal/"));

		map<TString, float> XSEC_ = {{"TT", 832.},
															  {"TTToHadronic", 377.96},
																{"TTToSemiLeptonic", 365.34},
																{"TTTo2L2Nu", 88.29}};

		XSECAll.insert(pair<TString, map<TString, float>>("2016", XSEC_));
		XSECAll.insert(pair<TString, map<TString, float>>("2017", XSEC_));
		XSECAll.insert(pair<TString, map<TString, float>>("2018", XSEC_));

		map<TString, TString> files2016 = {{"TTToHadronic", "2016/HistoCutFlowJetMassSoftDrop_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
		                                   {"TTToSemiLeptonic", "2016/HistoCutFlowJetMassSoftDrop_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
		                                   {"TTTo2L2Nu", "2016/HistoCutFlowJetMassSoftDrop_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root"}};

		map<TString, TString> files2017 = {{"TTToHadronic",  "2017/HistoCutFlowJetMassSoftDrop_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTToSemiLeptonic", "2017/HistoCutFlowJetMassSoftDrop_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTTo2L2Nu", "2017/HistoCutFlowJetMassSoftDrop_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

		map<TString, TString> files2018 = {{"TTToHadronic",  "2018/HistoCutFlowJetMassSoftDrop_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTToSemiLeptonic", "2018/HistoCutFlowJetMassSoftDrop_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTTo2L2Nu", "2018/HistoCutFlowJetMassSoftDrop_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

		files.insert(pair<TString, map<TString, TString>>("2016", files2016));
		files.insert(pair<TString, map<TString, TString>>("2017", files2017));
		files.insert(pair<TString, map<TString, TString>>("2018", files2018));
	}

	//these are btagging Working points for each year Medium WP
	floatConstants.insert(pair<TString, float>("btagWP2016", 0.6321));
	floatConstants.insert(pair<TString, float>("btagWP2017", 0.4941));
	floatConstants.insert(pair<TString, float>("btagWP2018", 0.4184));

	topTaggerConstants.insert(pair<TString, float>("topTagger2016", 0.2));
	topTaggerConstants.insert(pair<TString, float>("topTagger2017", 0.0));
	topTaggerConstants.insert(pair<TString, float>("topTagger2018", 0.1));

	luminosity.insert(pair<TString,float>("luminosity2016", 35920));
	luminosity.insert(pair<TString,float>("luminosity2017", 41530));
	luminosity.insert(pair<TString,float>("luminosity2018", 59740));

}
