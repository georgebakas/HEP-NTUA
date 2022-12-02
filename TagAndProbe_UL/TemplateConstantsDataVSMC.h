#include "TString.h"

map<TString, map<TString, TString>> mcFiles;
map<TString, map<TString, float>> XSECall;

map<TString, TString> dataFiles;


void initFilesXSEC()
{
//----------------------------------------------------------------------------------------------------------------
	//data files:
	dataFiles.insert(pair<TString, TString>("2016_preVFP","JetHT_Run2016-21Feb2020_UL2016_HIPM-v1.root"));
	dataFiles.insert(pair<TString, TString>("2016_postVFP","JetHT_Run2016-21Feb2020_UL2016-v1.root"));
	dataFiles.insert(pair<TString, TString>("2017","JetHT_Run2017-UL2017_MiniAODv2-v1.root"));
	dataFiles.insert(pair<TString, TString>("2018","JetHT_Run2018-UL2018_MiniAODv2-v1.root"));
//----------------------------------------------------------------------------------------------------------------	


//----------------------------------------------------------------------------------------------------------------
//ALL MC files together here:
	map<TString, TString> files_16pre = {
                                   {"QCD_HT300to500", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT500to700", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT700to1000", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1000to1500", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1500to2000", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT2000toInf", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
                                   {"TTToHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
								   {"TTToSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
								   {"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"},
                                   {"ST_tW_top_5f_NoFullyHadronicDecays", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
	                               {"ST_t-channel_top_5f_InclusiveDecays", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_t-channel_antitop_4f_InclusiveDecays", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
                                   {"ST_t-channel_top_4f_InclusiveDecays", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"}};
    
    map<TString, TString> files_16post = {
                                   {"QCD_HT300to500", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT500to700", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT700to1000", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1000to1500", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1500to2000", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT2000toInf", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
                                   {"TTToHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"TTToSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"},
                                   {"ST_tW_top_5f_NoFullyHadronicDecays", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
	                               {"ST_t-channel_top_5f_InclusiveDecays", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_t-channel_antitop_4f_InclusiveDecays", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
                                   {"ST_t-channel_top_4f_InclusiveDecays", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"}};

    map<TString, TString> files_17 = {
                                   {"QCD_HT300to500", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT500to700", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT700to1000", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1000to1500", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1500to2000", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT2000toInf", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
                                   {"TTToHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"TTToSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"},
                                   {"ST_tW_top_5f_NoFullyHadronicDecays", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_20UL.root"},
	                               {"ST_t-channel_top_5f_InclusiveDecays", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               //{"ST_t-channel_antitop_4f_InclusiveDecays", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
                                   {"ST_t-channel_top_4f_InclusiveDecays", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"}};

    map<TString, TString> files_18 = {
                                   {"QCD_HT300to500", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT500to700", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT700to1000", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1000to1500", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT1500to2000", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
	                               {"QCD_HT2000toInf", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL.root"},
                                   {"TTToHadronic", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"TTToSemiLeptonic", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"TTTo2L2Nu", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"},
                                   {"ST_tW_top_5f_NoFullyHadronicDecays", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_t-channel_top_5f_InclusiveDecays", "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_19UL.root"},
	                               {"ST_t-channel_antitop_4f_InclusiveDecays", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"},
                                   {"ST_t-channel_top_4f_InclusiveDecays", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_19UL.root"}};

	
    //qcd files MC:
    map<TString, float> XSEC16_pre =  {
                                   {"QCD_HT300to500",315400},
	                               {"QCD_HT500to700",32260},
	                               {"QCD_HT700to1000",6830},
	                               {"QCD_HT1000to1500",1207},
	                               {"QCD_HT1500to2000",119.1},
	                               {"QCD_HT2000toInf",25.16},
                                   {"TTToHadronic",377.96},
								   {"TTToSemiLeptonic",365.34},
								   {"TTTo2L2Nu",88.29},
                                   {"ST_tW_top_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_t-channel_top_5f_InclusiveDecays", 82.52},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", 119.7},
	                               {"ST_t-channel_antitop_4f_InclusiveDecays",35.6},
                                   {"ST_t-channel_top_4f_InclusiveDecays",35.6}};
    
    map<TString, float> XSEC16_post =  {{"QCD_HT300to500",351400},
	                               {"QCD_HT500to700",32260},
	                               {"QCD_HT700to1000",6830},
	                               {"QCD_HT1000to1500",1207},
	                               {"QCD_HT1500to2000",119.1},
	                               {"QCD_HT2000toInf",25.16},
                                   {"TTToHadronic",377.96},
								   {"TTToSemiLeptonic",365.34},
								   {"TTTo2L2Nu",88.29},
                                   {"ST_tW_top_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_t-channel_top_5f_InclusiveDecays", 82.52},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", 119.7},
	                               {"ST_t-channel_antitop_4f_InclusiveDecays",35.6},
                                   {"ST_t-channel_top_4f_InclusiveDecays",35.6}};

    map<TString, float> XSEC17 =    {{"QCD_HT300to500",322600},
	                               {"QCD_HT500to700",29980},
	                               {"QCD_HT700to1000",6334},
	                               {"QCD_HT1000to1500",1088},
	                               {"QCD_HT1500to2000",99.11},
	                               {"QCD_HT2000toInf",20.23},
                                   {"TTToHadronic",377.96},
								   {"TTToSemiLeptonic",365.34},
								   {"TTTo2L2Nu",88.29},
	                               {"ST_tW_top_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_t-channel_top_5f_InclusiveDecays", 82.52},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", 119.7},
	                               {"ST_t-channel_antitop_4f_InclusiveDecays",35.6},
                                   {"ST_t-channel_top_4f_InclusiveDecays",35.6}};

	map<TString, float> XSEC18 =    {{"QCD_HT300to500",323400},
	                               {"QCD_HT500to700",30140},
	                               {"QCD_HT700to1000",6310},
	                               {"QCD_HT1000to1500",1094},
	                               {"QCD_HT1500to2000",99.38},
	                               {"QCD_HT2000toInf",20.2},
                                   {"TTToHadronic",377.96},
								   {"TTToSemiLeptonic",365.34},
								   {"TTTo2L2Nu",88.29},
                                   {"ST_tW_top_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_tW_antitop_5f_NoFullyHadronicDecays", 38.09},
	                               {"ST_t-channel_top_5f_InclusiveDecays", 82.52},
	                               {"ST_t-channel_antitop_5f_InclusiveDecays", 119.7},
	                               {"ST_t-channel_antitop_4f_InclusiveDecays",35.6},
                                   {"ST_t-channel_top_4f_InclusiveDecays",35.6}};

	mcFiles.insert(pair<TString, map<TString, TString>>("2016_preVFP",files_16pre));
	mcFiles.insert(pair<TString, map<TString, TString>>("2016_postVFP",files_16post));
	mcFiles.insert(pair<TString, map<TString, TString>>("2017",files_17));
	mcFiles.insert(pair<TString, map<TString, TString>>("2018",files_18));

	XSECall.insert(pair<TString, map<TString, float>>("2016_preVFP",XSEC16_pre));
	XSECall.insert(pair<TString, map<TString, float>>("2016_postVFP",XSEC16_post));
	XSECall.insert(pair<TString, map<TString, float>>("2017",XSEC17));
	XSECall.insert(pair<TString, map<TString, float>>("2018",XSEC18));

//----------------------------------------------------------------------------------------------------------------

}
