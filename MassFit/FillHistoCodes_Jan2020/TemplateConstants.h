#include "TString.h"

map<TString, map<TString, TString>> files;
map<TString, map<TString, TString>> filesReduced;

map<TString, map<TString, TString>> mttFiles;
map<TString, map<TString, TString>> ttNominalFiles;
map<TString, map<TString, TString>> subdominantBkgFiles;
map<TString, map<TString, TString>> qcdBkgFiles;

map<TString, map<TString, float>> mttXSEC;
map<TString, map<TString, float>> ttNominalXSEC;
map<TString, map<TString, float>> subdominantBkgXSEC;
map<TString, map<TString, float>> qcdBkgXSEC;

map<TString, TString> dataFiles;
map<TString, float> Nbkg2Constants;
map<TString, float> Nbkg2ConstantsErrors;
map<TString, int> variableConstant;
TColor color;

map<TString, float> luminosity;
map<TString, float> luminosityCR;
map<TString, float> deepCSVFloatMap;
map<TString, float> topTaggerCuts;
map<TString, int> triggerSRConst;
map<TString, int> triggerCRConst;
TString eosDataPath, eosPathMC;

void initFilesMapping(bool isLoose)
{
//----------------------------------------------------------------------------------------------------------------
	//data files:
	eosDataPath = "/eos/cms/store/user/gbakas/ttbar/JetHT/ul-";
	dataFiles.insert(pair<TString, TString>("2016_preVFP","JetHT_Run2016-21Feb2020_UL2016_HIPM-v1.root"));
	dataFiles.insert(pair<TString, TString>("2016_postVFP","JetHT_Run2016-21Feb2020_UL2016-v1.root"));
	dataFiles.insert(pair<TString, TString>("2017","JetHT_Run2017-UL2017_MiniAODv2-v1.root"));
	dataFiles.insert(pair<TString, TString>("2018","JetHT_Run2018-UL2018_MiniAODv2-v1.root"));

//----------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------
	eosPathMC = "/eos/cms/store/user/gbakas/ttbar/topTagger/ul-";
	//qcd files MC:

	map<TString, TString> qcd16 = {{"300-500", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"500-700", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"700-1000", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"1000-1500", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"1500-2000", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"2000-Inf", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"}};

	map<TString, TString> qcd17 = {{"300-500", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"500-700", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"700-1000", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"1000-1500", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"1500-2000", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"2000-Inf", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"}};

	map<TString, TString> qcd18 = {{"300-500", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"500-700", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"700-1000", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"1000-1500", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"1500-2000", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"2000-Inf", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"}};

	map<TString, float> qcd16XSEC = {{"300-500",351400},
	                               {"500-700",32260},
	                               {"700-1000",6830},
	                               {"1000-1500",1207},
	                               {"1500-2000",119.1},
	                               {"2000-Inf",25.16}};

    map<TString, float> qcd17XSEC = {{"300-500",322600},
	                               {"500-700",29980},
	                               {"700-1000",6334},
	                               {"1000-1500",1088},
	                               {"1500-2000",99.11},
	                               {"2000-Inf",20.23}};

	map<TString, float> qcd18XSEC = {{"300-500",323400},
	                               {"500-700",30140},
	                               {"700-1000",6310},
	                               {"1000-1500",1094},
	                               {"1500-2000",99.38},
	                               {"2000-Inf",20.2}};

	qcdBkgFiles.insert(pair<TString, map<TString, TString>>("2016_preVFP", qcd16));
	qcdBkgFiles.insert(pair<TString, map<TString, TString>>("2016_postVFP", qcd16));
	qcdBkgFiles.insert(pair<TString, map<TString, TString>>("2017", qcd17));
	qcdBkgFiles.insert(pair<TString, map<TString, TString>>("2018", qcd18));

	qcdBkgXSEC.insert(pair<TString, map<TString, float>>("2016_preVFP",qcd16XSEC));
	qcdBkgXSEC.insert(pair<TString, map<TString, float>>("2016_postVFP",qcd16XSEC));
	qcdBkgXSEC.insert(pair<TString, map<TString, float>>("2017",qcd17XSEC));
	qcdBkgXSEC.insert(pair<TString, map<TString, float>>("2018",qcd18XSEC));
//----------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------
//subdominant bkgs:
	map<TString, TString> sub16_pre = {{"ST_tW_top_5f", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_tW_antitop_5f", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
									{"ST_t-channel_top_4f", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									{"ST_t-channel_antitop_4f", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									{"ST_t-channel_top_5f", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
									{"ST_t-channel_antitop_5f", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"WJetsToQQ_HT-200to400", "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-400to600", "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-600to800", "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-800toInf", "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root"}};

	map<TString, TString> sub16_post = {{"ST_tW_top_5f", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_tW_antitop_5f", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
									{"ST_t-channel_top_4f", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									{"ST_t-channel_antitop_4f", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									{"ST_t-channel_top_5f", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_t-channel_antitop_5f", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"WJetsToQQ_HT-200to400", "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-400to600", "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-600to800", "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-800toInf", "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root"}};

	map<TString, TString> sub17 = {{"ST_tW_top_5f", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_tW_antitop_5f", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
									{"ST_t-channel_top_4f", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									//{"ST_t-channel_antitop_4f", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									{"ST_t-channel_top_5f", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_t-channel_antitop_5f", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"WJetsToQQ_HT-200to400", "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-400to600", "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-600to800", "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-800toInf", "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root"}};

	map<TString, TString> sub18 = {{"ST_tW_top_5f", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_tW_antitop_5f", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_t-channel_top_4f", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									{"ST_t-channel_antitop_4f", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
									{"ST_t-channel_top_5f", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"ST_t-channel_antitop_5f", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
									{"WJetsToQQ_HT-200to400", "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-400to600", "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-600to800", "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
									{"WJetsToQQ_HT-800toInf", "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root"}};



	map<TString, float> sub16XSEC_pre = {{"ST_tW_top_5f", 38.09},
									{"ST_tW_antitop_5f", 38.09},
									{"ST_t-channel_top_4f", 35.6},
									{"ST_t-channel_antitop_4f", 35.6},
									{"ST_t-channel_top_5f", 119.7},
									{"ST_t-channel_antitop_5f", 82.52},
									{"WJetsToQQ_HT-200to400", 2549.0},
									{"WJetsToQQ_HT-400to600", 276.5},
									{"WJetsToQQ_HT-600to800", 59.25},
									{"WJetsToQQ_HT-800toInf", 28.75}};

	map<TString, float> sub16XSEC_post = {{"ST_tW_top_5f", 38.09},
									{"ST_tW_antitop_5f", 38.09},
									{"ST_t-channel_top_4f", 35.6},
									{"ST_t-channel_antitop_4f", 35.6},
									{"ST_t-channel_top_5f", 119.7},
									{"ST_t-channel_antitop_5f", 82.52},
									{"WJetsToQQ_HT-200to400", 2549.0},
									{"WJetsToQQ_HT-400to600", 276.5},
									{"WJetsToQQ_HT-600to800", 59.25},
									{"WJetsToQQ_HT-800toInf", 28.75}};

	map<TString, float> sub17XSEC = {{"ST_tW_top_5f", 34.91},
									{"ST_tW_antitop_5f", 34.97},
									{"ST_t-channel_top_4f", 115.3},
									//{"ST_t-channel_antitop_4f", 69.09},
									{"ST_t-channel_top_5f", 119.7},
									{"ST_t-channel_antitop_5f", 71.74},
									{"WJetsToQQ_HT-200to400", 2549.0},
									{"WJetsToQQ_HT-400to600", 276.5},
									{"WJetsToQQ_HT-600to800", 59.25},
									{"WJetsToQQ_HT-800toInf", 28.75}};

	map<TString, float> sub18XSEC = {{"ST_tW_top_5f", 34.91},
									{"ST_tW_antitop_5f", 34.97},
									{"ST_t-channel_top_4f", 115.3},
									{"ST_t-channel_antitop_4f", 69.09},
									{"ST_t-channel_top_5f", 119.7},
									{"ST_t-channel_antitop_5f", 71.74},
									{"WJetsToQQ_HT-200to400", 2549.0},
									{"WJetsToQQ_HT-400to600", 276.5},
									{"WJetsToQQ_HT-600to800", 59.25},
									{"WJetsToQQ_HT-800toInf", 28.75}};

	subdominantBkgFiles.insert(pair<TString, map<TString, TString>>("2016_preVFP",sub16_pre));
	subdominantBkgFiles.insert(pair<TString, map<TString, TString>>("2016_postVFP",sub16_post));
	subdominantBkgFiles.insert(pair<TString, map<TString, TString>>("2017",sub17));
	subdominantBkgFiles.insert(pair<TString, map<TString, TString>>("2018",sub18));

	subdominantBkgXSEC.insert(pair<TString, map<TString, float>>("2016_preVFP",sub16XSEC_pre));
	subdominantBkgXSEC.insert(pair<TString, map<TString, float>>("2016_postVFP",sub16XSEC_post));
	subdominantBkgXSEC.insert(pair<TString, map<TString, float>>("2017",sub17XSEC));
	subdominantBkgXSEC.insert(pair<TString, map<TString, float>>("2018",sub18XSEC));
//----------------------------------------------------------------------------------------------------------------
	//Nominal TT files (MC):
	map<TString, TString> eosNomTT16_preVFP = {{"TTHadronic","TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTSemiLeptonic","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTTo2L2Nu","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

	map<TString, TString> eosNomTT16_postVFP = {{"TTHadronic","TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTSemiLeptonic","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTTo2L2Nu","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

	map<TString, TString> eosNomTT17 = {{"TTHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

	map<TString, TString> eosNomTT18 = {{"TTHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
									{"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};

	map<TString, float> eosNomTTXSEC16 = {{"TTHadronic",377.96},
									{"TTSemiLeptonic",365.34},
									{"TTTo2L2Nu",88.29}};

    map<TString, float> eosNomTTXSEC17 = {{"TTHadronic",377.96},
									{"TTSemiLeptonic",365.34},
									{"TTTo2L2Nu",88.29}};

    map<TString, float> eosNomTTXSEC18 = {{"TTHadronic",377.96},
									{"TTSemiLeptonic",365.34},
									{"TTTo2L2Nu",88.29}};

    ttNominalFiles.insert(pair<TString, map<TString, TString>>("2016_preVFP", eosNomTT16_preVFP));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2016_postVFP", eosNomTT16_postVFP));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2017", eosNomTT17));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2018", eosNomTT18));

	ttNominalXSEC.insert(pair<TString, map<TString, float>>("2016_preVFP",eosNomTTXSEC16));
	ttNominalXSEC.insert(pair<TString, map<TString, float>>("2016_postVFP",eosNomTTXSEC16));
	ttNominalXSEC.insert(pair<TString, map<TString, float>>("2017",eosNomTTXSEC17));
	ttNominalXSEC.insert(pair<TString, map<TString, float>>("2018",eosNomTTXSEC18));

//----------------------------------------------------------------------------------------------------------------


	//thsese will be used either for naming the files for FillHistograms, or either for getting the histo names
	map<TString, TString> files2016_pre = {{"data",  "2016_preVFP/Histo_Data_2016_preVFP_100.root"},
	                                       {"mcSig", "2016_preVFP/Histo_TT_Mtt-700toInf_100.root"},
	                                       {"mcSub", "2016_preVFP/Histo_SubdominantBkgs_100.root"},
	                               	       {"qcd"  , "2016_preVFP/Histo_QCD_HT300toInf_100.root"}};

	map<TString, TString> files2016Loose_pre = {{"data",  "2016_preVFP/Histo_Data_2016_preVFP_100_Loose.root"},
	                                   		    {"mcSig", "2016_preVFP/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   		    {"mcSub", "2016_preVFP/Histo_SubdominantBkgs_100_Loose.root"},
	                                   			{"qcd"  , "2016_preVFP/Histo_QCD_HT300toInf_100_Loose.root"}};

	map<TString, TString> files2016_post = {{"data",  "2016_postVFP/Histo_Data_2016_preVFP_100.root"},
	                                   		{"mcSig", "2016_postVFP/Histo_TT_Mtt-700toInf_100.root"},
	                                   		{"mcSub", "2016_postVFP/Histo_SubdominantBkgs_100.root"},
	                               	  		{"qcd"  , "2016_postVFP/Histo_QCD_HT300toInf_100.root"}};

	map<TString, TString> files2016Loose_post = {{"data",  "2016_postVFP/Histo_Data_2016_postVFP_100_Loose.root"},
	                                   			 {"mcSig", "2016_postVFP/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   			 {"mcSub", "2016_postVFP/Histo_SubdominantBkgs_100_Loose.root"},
	                                   			 {"qcd"  , "2016_postVFP/Histo_QCD_HT300toInf_100_Loose.root"}};

	map<TString, TString> files2017 = {{"data",  "2017/Histo_Data_2017_100.root"},
	                                   {"mcSig", "2017/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2017/Histo_SubdominantBkgs_100.root"},
	                                   {"qcd"  , "2017/Histo_QCD_HT300toInf_100.root"}};


    map<TString, TString> files2017Loose = {{"data",  "2017/Histo_Data_2017_100_Loose.root"},
	                                   		{"mcSig", "2017/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   		{"mcSub", "2017/Histo_SubdominantBkgs_100_Loose.root"},
	                                   		{"qcd"  , "2017/Histo_QCD_HT300toInf_100_Loose.root"}};

	map<TString, TString> files2018 = {{"data",  "2018/Histo_Data_2018_100.root"},
	                                   {"mcSig", "2018/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2018/Histo_SubdominantBkgs_100.root"},
	                                   {"qcd"  , "2018/Histo_QCD_HT300toInf_100.root"}};

	map<TString, TString> files2018Loose = {{"data",  "2018/Histo_Data_2018_100_Loose.root"},
	                                   		{"mcSig", "2018/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   		{"mcSub", "2018/Histo_SubdominantBkgs_100_Loose.root"},
	                                   		{"qcd"  , "2018/Histo_QCD_HT300toInf_100_Loose.root"}};

	files.insert(pair<TString, map<TString, TString>>("2016_preVFP", files2016_pre));
	files.insert(pair<TString, map<TString, TString>>("2016_postVFP", files2016_post));
	files.insert(pair<TString, map<TString, TString>>("2017", files2017));
	files.insert(pair<TString, map<TString, TString>>("2018", files2018));
	files.insert(pair<TString, map<TString, TString>>("2016_preVFP_Loose", files2016Loose_pre));
	files.insert(pair<TString, map<TString, TString>>("2016_postVFP_Loose", files2016Loose_post));
	files.insert(pair<TString, map<TString, TString>>("2017_Loose", files2017Loose));
	files.insert(pair<TString, map<TString, TString>>("2018_Loose", files2018Loose));

	map<TString, TString> files2016Reduced_pre = {{"data", "2016_preVFP/Histo_Data_2016_preVFP_100_reduced.root"},
	                                   		  {"mcSig", "2016_preVFP/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016_preVFP/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2016_preVFP/Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2016Reduced_post = {{"data", "2016_postVFP/Histo_Data_2016_postVFP_100_reduced.root"},
	                                   		  {"mcSig", "2016_postVFP/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016_postVFP/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2016_postVFP/Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2017Reduced = {{"data",  "2017/Histo_Data_2017_100_reduced.root"},
	                                   		  {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "20167Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2018Reduced = {{"data",  "2018/Histo_Data_2018_100_reduced.root"},
	                                   	      {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2018/Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2016ReducedLoose_pre = {{"data", "2016_preVFP/Histo_Data_2016_preVFP_100_reduced_Loose.root"},
	                                   		       {"mcSig", "2016_preVFP/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		       {"mcSub", "2016_preVFP/Histo_SubdominantBkgs_100_reduced_Loose.root"},
	                                   		       {"qcd"  , "2016_preVFP/Histo_QCD_HT300toInf_100_reduced_Loose.root"}};

	map<TString, TString> files2016ReducedLoose_post = {{"data", "2016_postVFP/Histo_Data_2016_postVFP_100_reduced_Loose.root"},
	                                   		       {"mcSig", "2016_postVFP/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		       {"mcSub", "2016_postVFP/Histo_SubdominantBkgs_100_reduced_Loose.root"},
	                                   		       {"qcd"  , "2016_postVFP/Histo_QCD_HT300toInf_100_reduced_Loose.root"}};

	map<TString, TString> files2017ReducedLoose = {{"data",  "2017/Histo_Data_2017_100_reduced_Loose.root"},
	                                   		  	   {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		  	   {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced_Loose.root"},
	                                   		  	   {"qcd"  , "2017/Histo_QCD_HT300toInf_100_reduced_Loose.root"}};

	map<TString, TString> files2018ReducedLoose = {{"data",  "2018/Histo_Data_2018_100_reduced_Loose.root"},
	                                   	           {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		       {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced_Loose.root"},
	                                   		       {"qcd"  , "2018/Histo_QCD_HT300toInf_100_reduced_Loose.root"}};

	filesReduced.insert(pair<TString, map<TString, TString>>("2016_preVFP", files2016Reduced_pre));
	filesReduced.insert(pair<TString, map<TString, TString>>("2016_postVFP", files2016Reduced_post));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017", files2017Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018", files2018Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2016_preVFP_Loose", files2016ReducedLoose_pre));
	filesReduced.insert(pair<TString, map<TString, TString>>("2016_postVFP_Loose", files2016ReducedLoose_post));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017_Loose", files2017ReducedLoose));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018_Loose", files2018ReducedLoose));

//---------------------------------------------------------------------------------------------------------------------------------

	luminosity.insert(pair<TString, float>("2016_preVFP", 19500));
	luminosity.insert(pair<TString, float>("2016_postVFP", 16500));
	luminosity.insert(pair<TString, float>("2017", 41480));
	luminosity.insert(pair<TString, float>("2018", 59830));

	luminosityCR.insert(pair<TString, float>("2016_preVFP", 1134));
	//luminosityCR.insert(pair<TString, float>("2016_preVFP", 904.58));
	luminosityCR.insert(pair<TString, float>("2016_postVFP", 564));
	//luminosityCR.insert(pair<TString, float>("2016_postVFP", 765.42));
	luminosityCR.insert(pair<TString, float>("2017", 41480));
	luminosityCR.insert(pair<TString, float>("2018", 59830));

	topTaggerCuts.insert(pair<TString, float>("2016_preVFP", 0.2));
	topTaggerCuts.insert(pair<TString, float>("2016_postVFP", 0.2));
	topTaggerCuts.insert(pair<TString, float>("2017", 0.0));
	topTaggerCuts.insert(pair<TString, float>("2018", 0.1));

	if(!isLoose)
	{
		deepCSVFloatMap.insert(pair<TString, float>("2016_preVFP", 0.6001));
		deepCSVFloatMap.insert(pair<TString, float>("2016_postVFP", 0.5847));
		deepCSVFloatMap.insert(pair<TString, float>("2017",0.4506));
		deepCSVFloatMap.insert(pair<TString, float>("2018",0.4168));
	}
	else
	{
		deepCSVFloatMap.insert(pair<TString, float>("2016_preVFP",0.2217));
		deepCSVFloatMap.insert(pair<TString, float>("2016_postVFP",0.2217));
		deepCSVFloatMap.insert(pair<TString, float>("2017",0.1355));
		deepCSVFloatMap.insert(pair<TString, float>("2018",0.1208));
	}

	triggerSRConst.insert(pair<TString, int>("2016_preVFP",2));
	triggerSRConst.insert(pair<TString, int>("2016_postVFP",2));
	triggerSRConst.insert(pair<TString, int>("2017",5));
	triggerSRConst.insert(pair<TString, int>("2018",5));

	triggerCRConst.insert(pair<TString, int>("2016_preVFP",4));
	triggerCRConst.insert(pair<TString, int>("2016_postVFP",4));
	triggerCRConst.insert(pair<TString, int>("2017",5));
	triggerCRConst.insert(pair<TString, int>("2018",5));

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
