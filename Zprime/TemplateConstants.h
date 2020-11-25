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


	map<TString, float> XSEC2016 = {{"M-1000_W-10", 3.205},
																	{"M-1000_W-100", 0.314},
																	{"M-1000_W-300", 0.09429},

																	{"M-2000_W-20", 0.1451},
																	{"M-2000_W-200", 0.01592},
																	{"M-2000_W-600", 0.005714},

																	{"M-2500_W-25", 0.04169},
																	{"M-2500_W-250", 0.005061},

																	{"M-3000_W-30", 0.01309},
																	{"M-3000_W-300", 0.001809},
																	{"M-3000_W-900", 0.0008196},

																	{"M-4000_W-40", 0.001568},
																	{"M-4000_W-400", 0.0003247},
																	{"M-4000_W-1200", 0.000193},

																	{"M-5000_W-50", 0.0002212},
																	{"M-5000_W-500", 8.866e-05},
																	{"M-5000_W-1500", 6.499e-05}};

	map<TString, float> XSEC2017 = {{"M1000_W10", 3.515},
																	{"M1000_W100", 0.3391},
																	{"M1000_W300", 0.1045},

																	{"M2000_W20", 0.1617},
																	{"M2000_W200", 0.01842},
																	{"M2000_W600", 0.006504},

																	{"M2500_W25", 0.04753},
																	{"M2500_W250", 0.00563},
																	{"M2500_W750", 0.002282},

																	{"M3000_W-30", 0.01505},
																	{"M3000_W-300", 0.002006},
																	{"M3000_W-900", 0.0009163},

																	{"M3500_W-35", 0.005014},
																	{"M3500_W-350", 0.0008262},
																	{"M3500_W-1050", 0.0004251},

																	{"M4000_W40", 0.001917},
																	{"M4000_W400", 0.0003791},
																	{"M4000_W1200", 0.0002185},

																	{"M5000_W50", 0.0003224},
																	{"M5000_W500", 9.517e-05},
																	{"M5000_W1500", 7.394e-05}};

	map<TString, float> XSEC2018 = {{"M1000_W10", 3.52},
																	{"M1000_W100", 0.3463},
																	{"M1000_W300", 0.1036},

																	{"M2000_W20", 0.165},
																	{"M2000_W200", 0.01825},
																	{"M2000_W600", 0.006399},

																	{"M2500_W25", 0.04715},
																	{"M2500_W250", 0.005708},
																	{"M2500_W750", 0.002255},

																	{"M3000_W30", 0.01495},
																	{"M3000_W300", 0.002056},
																	{"M3000_W900", 0.0009167},

																	{"M3500_W-35", 0.005108},
																	{"M3500_W-350", 0.0008352},
																	{"M3500_W-1050", 0.0004248},

																	{"M4000_W40", 0.001908},
																	{"M4000_W400", 0.0003779},
																	{"M4000_W1200", 0.0002186},

																	{"M5000_W50", 0.0003225},
																	{"M5000_W500", 0.0001055},
																	{"M5000_W1500", 7.46e-05}};


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
	//floatConstants.insert(pair<TString, float>("btagWP2016", 0.6321));
	//floatConstants.insert(pair<TString, float>("btagWP2017", 0.4941));
	//floatConstants.insert(pair<TString, float>("btagWP2018", 0.4184));

	//loose WP
	floatConstants.insert(pair<TString, float>("btagWP2016", 0.2217));
	floatConstants.insert(pair<TString, float>("btagWP2017", 0.1522));
	floatConstants.insert(pair<TString, float>("btagWP2018", 0.1241));

	//ttbar analysis WP
	//topTaggerConstants.insert(pair<TString, float>("topTagger2016", 0.2));
	//topTaggerConstants.insert(pair<TString, float>("topTagger2017", 0.0));
	//topTaggerConstants.insert(pair<TString, float>("topTagger2018", 0.1));

	//Z prime analysis WP
	topTaggerConstants.insert(pair<TString, float>("topTagger2016", 0.1));
	topTaggerConstants.insert(pair<TString, float>("topTagger2017", -0.1));
	topTaggerConstants.insert(pair<TString, float>("topTagger2018", 0.0));

	luminosity.insert(pair<TString,float>("luminosity2016", 35920));
	luminosity.insert(pair<TString,float>("luminosity2017", 41530));
	luminosity.insert(pair<TString,float>("luminosity2018", 59740));

}
