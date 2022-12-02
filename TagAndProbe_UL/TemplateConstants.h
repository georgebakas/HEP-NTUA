#include "TString.h"

map<TString, map<TString, float>> XSECAll;

map<TString, float> floatConstants;
map<TString, float> topTaggerConstants;
map<TString, float> luminosity;
map<TString, float> luminosityCR;
map<TString, TString> eospath;
map<TString, TString> eospath_nom;
map<int, TString> ps_weights;
map<int, TString> pdf_weights;
map<int, TString> scale_weights;
map<TString, int> variableConstant;
map<int, TString> vars;
map<TString, map<TString, TString>> ttNominalFiles;
map<TString, float> ttbarSigStrength;

void initFilesMapping()
{

	eospath.insert(pair<TString,TString>("2016_preVFP","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_preVFP/variations/"));
	eospath.insert(pair<TString,TString>("2016_postVFP","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_postVFP/variations/"));
	eospath.insert(pair<TString,TString>("2017","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/variations/"));
	eospath.insert(pair<TString,TString>("2018","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2018/variations/"));

	eospath_nom.insert(pair<TString,TString>("2016_preVFP","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_preVFP/Signal/"));
	eospath_nom.insert(pair<TString,TString>("2016_postVFP","/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2016_postVFP/Signal/"));
	eospath_nom.insert(pair<TString,TString>("2017","/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2017/"));
	eospath_nom.insert(pair<TString,TString>("2018","/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2018/"));

	map<TString, float> XSEC_ = {{"TT", 832.},
								{"TTJets", 832.},
								{"TTToHadronic", 377.96},
								{"TTToSemiLeptonic", 365.34},
								{"TTTo2L2Nu", 88.29}};

	XSECAll.insert(pair<TString, map<TString, float>>("2016_preVFP", XSEC_));
	XSECAll.insert(pair<TString, map<TString, float>>("2016_postVFP", XSEC_));
	XSECAll.insert(pair<TString, map<TString, float>>("2017", XSEC_));
	XSECAll.insert(pair<TString, map<TString, float>>("2018", XSEC_));

	//these are btagging Working points for each year Medium WP
	floatConstants.insert(pair<TString, float>("btagWP2016_preVFP", 0.6001));
	floatConstants.insert(pair<TString, float>("btagWP2016_postVFP", 0.5847));
	floatConstants.insert(pair<TString, float>("btagWP2017",0.4506));
	floatConstants.insert(pair<TString, float>("btagWP2018",0.4168));

	topTaggerConstants.insert(pair<TString, float>("topTagger2016_preVFP", 0.2));
	topTaggerConstants.insert(pair<TString, float>("topTagger2016_postVFP", 0.2));
	topTaggerConstants.insert(pair<TString, float>("topTagger2017", 0.0));
	topTaggerConstants.insert(pair<TString, float>("topTagger2018", 0.1));

	luminosity.insert(pair<TString, float>("luminosity2016_preVFP", 19500));
	luminosity.insert(pair<TString, float>("luminosity2016_postVFP", 16500));
	luminosity.insert(pair<TString, float>("luminosity2017", 41480));
	luminosity.insert(pair<TString, float>("luminosity2018", 59830));
	luminosity.insert(pair<TString, float>("luminosityAll", 137310));


	luminosityCR.insert(pair<TString, float>("luminosity2016_preVFP", 1134));
	luminosityCR.insert(pair<TString, float>("luminosity2016_postVFP", 564));
	luminosityCR.insert(pair<TString, float>("luminosity2017", 41480));
	luminosityCR.insert(pair<TString, float>("luminosity2018", 59830));

	//ps weights
	ps_weights.insert(pair<int, TString>(0, "nom0"));
	ps_weights.insert(pair<int, TString>(1, "nom1"));
	ps_weights.insert(pair<int, TString>(2, "isrRedHi"));
	ps_weights.insert(pair<int, TString>(3, "fsrRedHi"));
	ps_weights.insert(pair<int, TString>(4, "isrRedLo"));
	ps_weights.insert(pair<int, TString>(5, "fsrRedLo"));
	ps_weights.insert(pair<int, TString>(6, "isrDefHi"));
	ps_weights.insert(pair<int, TString>(7, "fsrDefHi"));
	ps_weights.insert(pair<int, TString>(8, "isrDefLo"));
	ps_weights.insert(pair<int, TString>(9, "fsrDefLo"));
	ps_weights.insert(pair<int, TString>(10, "isrConHi"));
	ps_weights.insert(pair<int, TString>(11, "fsrConHi"));
	ps_weights.insert(pair<int, TString>(12, "isrConLo"));
	ps_weights.insert(pair<int, TString>(13, "fsrConLo"));

	for(int i =0; i<101; i++)
	{
		pdf_weights.insert(pair<int, TString>(i, TString::Format("pdf_%d",i)));
	}

	scale_weights.insert(pair<int, TString>(0, "scale_2"));
	scale_weights.insert(pair<int, TString>(1, "scale_3"));
	scale_weights.insert(pair<int, TString>(2, "scale_4"));
	scale_weights.insert(pair<int, TString>(3, "scale_5"));
	scale_weights.insert(pair<int, TString>(4, "scale_7"));
	scale_weights.insert(pair<int, TString>(5, "scale_9"));

	variableConstant.insert(pair<TString, int>("mJJ",  0));
	variableConstant.insert(pair<TString, int>("ptJJ", 1));
	variableConstant.insert(pair<TString, int>("yJJ",  2));
	variableConstant.insert(pair<TString, int>("jetPt0", 3));
	variableConstant.insert(pair<TString, int>("jetPt1", 4));
	variableConstant.insert(pair<TString, int>("jetY0", 5));
	variableConstant.insert(pair<TString, int>("jetY1", 6));
	variableConstant.insert(pair<TString, int>("chi", 7));
	variableConstant.insert(pair<TString, int>("cosTheta_0", 8));
	variableConstant.insert(pair<TString, int>("cosTheta_1", 9));
	

	vars.insert(pair<int, TString>(0, "mJJ"));
	vars.insert(pair<int, TString>(1, "ptJJ"));
	vars.insert(pair<int, TString>(2, "yJJ"));
	vars.insert(pair<int, TString>(3, "jetPt0"));
	vars.insert(pair<int, TString>(4, "jetPt1"));
	vars.insert(pair<int, TString>(5, "jetY0"));
	vars.insert(pair<int, TString>(6, "jetY1"));
	vars.insert(pair<int, TString>(7, "chi"));
	vars.insert(pair<int, TString>(8, "cosTheta_0"));
	vars.insert(pair<int, TString>(9, "cosTheta_1"));
	

	// Only for JES

	map<TString, TString> eosNomTT16 = {{"TTToHadronic","TTToHadronic_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
									   	{"TTToSemiLeptonic","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
									   	{"TTTo2L2Nu","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_20UL.root"}};
	/*map<TString, TString> eosNomTT16_post = {{"Hadronic","TTToHadronic_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
																			{"SemiLeptonic","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
																			{"Dilepton","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_20UL.root"}};*/

	map<TString, TString> eosNomTT17 = {{"TTToHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               	   {"TTToSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               	   {"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_19UL.root"}};

	/* map<TString, TString> eosNomTT18 = {{"TTToHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
	                               	   {"TTToSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
	                               	   {"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_20UL.root"}}; */

	map<TString, TString> eosNomTT18 = {{"TTToHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	   {"TTToSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	   {"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2016_preVFP", eosNomTT16));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2016_postVFP", eosNomTT16));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2017", eosNomTT17));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2018", eosNomTT18));

	ttbarSigStrength.insert(pair<TString,float>("2016_preVFP", 0.691));
	ttbarSigStrength.insert(pair<TString,float>("2016_postVFP", 0.640));
	ttbarSigStrength.insert(pair<TString,float>("2017", 0.665));
	ttbarSigStrength.insert(pair<TString,float>("2018", 0.675));



}
