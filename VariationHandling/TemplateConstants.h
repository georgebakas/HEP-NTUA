#include "TString.h"

map<TString, map<TString, float>> XSECAll;

map<TString, float> floatConstants;
map<TString, float> topTaggerConstants;
map<TString, float> luminosity;
map<TString, TString> eospath;

void initFilesMapping()
{

	eospath.insert(pair<TString,TString>("2016","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/variations/combined/"));
	eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/variations/combined/"));
	eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/variations/combined/"));

	map<TString, float> XSEC_ = {{"TT", 832.},
														  {"TTToHadronic", 377.96},
															{"TTToSemiLeptonic", 365.34},
															{"TTTo2L2Nu", 88.29}};

	XSECAll.insert(pair<TString, map<TString, float>>("2016", XSEC_));
	XSECAll.insert(pair<TString, map<TString, float>>("2017", XSEC_));
	XSECAll.insert(pair<TString, map<TString, float>>("2018", XSEC_));

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
