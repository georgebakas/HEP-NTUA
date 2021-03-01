#include "TString.h"

map<TString, map<TString, TString>> files;
map<TString, map<TString, float>> XSECAll;

map<TString, float> floatConstants;
map<TString, float> topTaggerConstants;
map<TString, float> luminosity;
map<TString, TString> eospath;

void initFilesMapping(bool isTTbar)
{
  if(!isTTbar)
	{
	eospath.insert(pair<TString,TString>("2016_preVFP","/eos/cms/store/user/gbakas/ZprimeToTT/ul-2016_preVFP/"));
	eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ZprimeToTT/ul-2017/"));
	eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ZprimeToTT/ul-2018/"));


	map<TString, float> XSEC2016 = {{"M1000_W10", 3.205},
									{"M1000_W100", 0.314},
									{"M1000_W300", 0.09429},

									{"M2000_W20", 0.1451},
									{"M2000_W200", 0.01592},
									{"M2000_W600", 0.005714},

									{"M2500_W25", 0.04169},
									{"M2500_W250", 0.005061},

									{"M3000_W30", 0.01309},
									{"M3000_W300", 0.001809},
									{"M3000_W900", 0.0008196},

									{"M3500_W35", 0.005014},
									{"M3500_W350", 0.0008262},
									{"M3500_W1050", 0.0004251},

									{"M4000_W40", 0.001568},
									{"M4000_W400", 0.0003247},
									{"M4000_W1200", 0.000193},

									{"M5000_W50", 0.0002212},
									{"M5000_W500", 8.866e-05},
									{"M5000_W1500", 6.499e-05}};

	map<TString, float> XSEC2017 = {{"M1000_W10", 3.515},
									{"M1000_W100", 0.3391},
									{"M1000_W300", 0.1045},

									{"M2000_W20", 0.1617},
									{"M2000_W200", 0.01842},
									{"M2000_W600", 0.006504},

									{"M2500_W25", 0.04753},
									{"M2500_W250", 0.00563},
									{"M2500_W750", 0.002282},

									{"M3000_W30", 0.01505},
									{"M3000_W300", 0.002006},
									{"M3000_W900", 0.0009163},

									{"M3500_W35", 0.005014},
									{"M3500_W350", 0.0008262},
									{"M3500_W1050", 0.0004251},

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

									{"M3500_W35", 0.005108},
									{"M3500_W350", 0.0008352},
									{"M3500_W1050", 0.0004248},

									{"M4000_W40", 0.001908},
									{"M4000_W400", 0.0003779},
									{"M4000_W1200", 0.0002186},

									{"M5000_W50", 0.0003225},
									{"M5000_W500", 0.0001055},
									{"M5000_W1500", 7.46e-05}};


 	XSECAll.insert(pair<TString, map<TString, float>>("2016_preVFP", XSEC2016));
	XSECAll.insert(pair<TString, map<TString, float>>("2016_postVFP", XSEC2016));
	XSECAll.insert(pair<TString, map<TString, float>>("2017", XSEC2017));
	XSECAll.insert(pair<TString, map<TString, float>>("2018", XSEC2018));

	}
	else{
		eospath.insert(pair<TString,TString>("2016_preVFP","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_preVFP/Signal/"));
		eospath.insert(pair<TString,TString>("2016_postVFP","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_postVFP/Signal/"));
		eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Signal/"));
		eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2018/Signal/"));

		map<TString, float> XSEC_ = {{"TT", 832.},
									{"TTToHadronic", 377.96},
									{"TTToSemiLeptonic", 365.34},
									{"TTTo2L2Nu", 88.29}};

		XSECAll.insert(pair<TString, map<TString, float>>("2016_preVFP", XSEC_));
		XSECAll.insert(pair<TString, map<TString, float>>("2016_postVFP", XSEC_));
		XSECAll.insert(pair<TString, map<TString, float>>("2017", XSEC_));
		XSECAll.insert(pair<TString, map<TString, float>>("2018", XSEC_));

		map<TString, TString> files2016_pre = {{"TTToHadronic", "2016_preVFP/HistoCutFlowJetMassSoftDrop_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
		                                   {"TTToSemiLeptonic", "2016_preVFP/HistoCutFlowJetMassSoftDrop_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
		                                   {"TTTo2L2Nu", "2016_preVFP/HistoCutFlowJetMassSoftDrop_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root"}};

		map<TString, TString> files2016_post = {{"TTToHadronic", "2016_postVFP/HistoCutFlowJetMassSoftDrop_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
		                                   {"TTToSemiLeptonic", "2016_postVFP/HistoCutFlowJetMassSoftDrop_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
		                                   {"TTTo2L2Nu", "2016_postVFP/HistoCutFlowJetMassSoftDrop_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root"}};


		map<TString, TString> files2017 = {{"TTToHadronic",  "2017/HistoCutFlowJetMassSoftDrop_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTToSemiLeptonic", "2017/HistoCutFlowJetMassSoftDrop_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTTo2L2Nu", "2017/HistoCutFlowJetMassSoftDrop_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

		map<TString, TString> files2018 = {{"TTToHadronic",  "2018/HistoCutFlowJetMassSoftDrop_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTToSemiLeptonic", "2018/HistoCutFlowJetMassSoftDrop_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
		                                   {"TTTo2L2Nu", "2018/HistoCutFlowJetMassSoftDrop_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

		files.insert(pair<TString, map<TString, TString>>("2016_preVFP", files2016_pre));
		files.insert(pair<TString, map<TString, TString>>("2016_postVFP", files2016_post));
		files.insert(pair<TString, map<TString, TString>>("2017", files2017));
		files.insert(pair<TString, map<TString, TString>>("2018", files2018));
	}

	//these are btagging Working points for each year Medium WP
	floatConstants.insert(pair<TString, float>("btagWP2016_preVFP", 0.6321));
	floatConstants.insert(pair<TString, float>("btagWP2016_postVFP", 0.6321));
	floatConstants.insert(pair<TString, float>("btagWP2017", 0.4941));
	floatConstants.insert(pair<TString, float>("btagWP2018", 0.4184));

	//ttbar analysis WP
	topTaggerConstants.insert(pair<TString, float>("topTagger2016_preVFP", 0.2));
	topTaggerConstants.insert(pair<TString, float>("topTagger2016_postVFP", 0.2));
	topTaggerConstants.insert(pair<TString, float>("topTagger2017", 0.0));
	topTaggerConstants.insert(pair<TString, float>("topTagger2018", 0.1));

	luminosity.insert(pair<TString,float>("luminosity2016_preVFP", 19500));
	luminosity.insert(pair<TString,float>("luminosity2016_postVFP", 16500));
	luminosity.insert(pair<TString,float>("luminosity2017", 41480));
	luminosity.insert(pair<TString,float>("luminosity2018", 59830));

}
