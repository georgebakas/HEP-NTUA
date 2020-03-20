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
	eosPathMC = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-";
	//files to be used from eos:
//----------------------------------------------------------------------------------------------------------------	
	//Mtt files (MC):
	map<TString, TString> eosMtt16 = {{"700-1000", "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"},
	                                  {"1000-Inf", "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"}};
	
    map<TString, TString> eosMtt17 = {{"700-1000", "TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8.root"},
	                                  {"1000-Inf", "TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root"}};	                              

    map<TString, TString> eosMtt18 = {{"700-1000", "TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8.root"},
	                                  {"1000-Inf", "TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root"}};	
	
	map<TString, float> eosMttXSEC16 = {{"700-1000",69.64},
	                                  {"1000-Inf",16.74}};
	
    map<TString, float> eosMttXSEC17 = {{"700-1000",69.64},
	                                  {"1000-Inf",16.74}};	                              

    map<TString, float> eosMttXSEC18 = {{"700-1000",69.64},
	                                  {"1000-Inf",16.74}};		                                  

    mttFiles.insert(pair<TString, map<TString, TString>>("2016", eosMtt16));
	mttFiles.insert(pair<TString, map<TString, TString>>("2017", eosMtt17));
	mttFiles.insert(pair<TString, map<TString, TString>>("2018", eosMtt18));                               

	mttXSEC.insert(pair<TString, map<TString, float>>("2016",eosMttXSEC16));
	mttXSEC.insert(pair<TString, map<TString, float>>("2017",eosMttXSEC17));
	mttXSEC.insert(pair<TString, map<TString, float>>("2018",eosMttXSEC18));
//----------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------
	//data files:
	eosDataPath = "/eos/cms/store/user/gbakas/ttbar/JetHT/";
	dataFiles.insert(pair<TString, TString>("2016","JetHT_Run2016-17Jul2018.root"));
	dataFiles.insert(pair<TString, TString>("2017","JetHT_Run2017-31Mar2018-v1.root"));
	dataFiles.insert(pair<TString, TString>("2018","JetHT_Run2018-17Sep2018-v1.root"));

//----------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------
	//qcd files MC:
	map<TString, TString> qcd16 = {{"300-500", "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                               {"500-700", "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                               {"700-1000", "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                               {"1000-1500", "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                               {"1500-2000", "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"},
	                               {"2000-Inf", "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"}};
	
	map<TString, TString> qcd17 = {{"300-500", "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root"},
	                               {"500-700", "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.root"},
	                               {"700-1000", "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.root"},
	                               {"1000-1500", "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.root"},
	                               {"1500-2000", "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root"},
	                               {"2000-Inf", "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.root"}};

	map<TString, TString> qcd18 = {{"300-500", "QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"500-700", "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"700-1000", "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"1000-1500", "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"1500-2000", "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"2000-Inf", "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root"}};

	map<TString, float> qcd16XSEC =    {{"300-500",351400},
	                               {"500-700",32260},
	                               {"700-1000",6830},
	                               {"1000-1500",1207},
	                               {"1500-2000",119.1},
	                               {"2000-Inf",25.16}};	

    map<TString, float> qcd17XSEC =    {{"300-500",322600},
	                               {"500-700",29980},
	                               {"700-1000",6334},
	                               {"1000-1500",1088},
	                               {"1500-2000",99.11},
	                               {"2000-Inf",20.23}};	
	
	map<TString, float> qcd18XSEC =    {{"300-500",323400},
	                               {"500-700",30140},
	                               {"700-1000",6310},
	                               {"1000-1500",1094},
	                               {"1500-2000",99.38},
	                               {"2000-Inf",20.2}};		                               

	qcdBkgFiles.insert(pair<TString, map<TString, TString>>("2016", qcd16));
	qcdBkgFiles.insert(pair<TString, map<TString, TString>>("2017", qcd17));
	qcdBkgFiles.insert(pair<TString, map<TString, TString>>("2018", qcd18)); 

	qcdBkgXSEC.insert(pair<TString, map<TString, float>>("2016",qcd16XSEC));
	qcdBkgXSEC.insert(pair<TString, map<TString, float>>("2017",qcd17XSEC));
	qcdBkgXSEC.insert(pair<TString, map<TString, float>>("2018",qcd18XSEC));
//----------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------
//subdominant bkgs:
	map<TString, TString> sub16 = {{"DY", "DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root"},
	                               {"WJetsToQQ_HT180", "WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root"},
	                               {"ST_tW_top_5f", "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root"},
	                               {"ST_tW_antitop_5f", "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root"},
	                               {"ST_t-channel_top_4f", "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root"},
	                               {"ST_t-channel_antitop_4f", "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root"}};

    map<TString, TString> sub17 = {{"DY", "DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8.root"},
	                               {"WJetsToQQ_HT400", "WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"WJetsToQQ_HT600", "WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"ST_tW_top_5f", "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"ST_tW_antitop_5f", "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"ST_t-channel_top_4f", "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root"},
	                               {"ST_t-channel_antitop_4f", "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root"}};	   

    map<TString, TString> sub18 = {{"DY", "DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8.root"},
	                               {"WJetsToQQ_HT400", "WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"WJetsToQQ_HT600", "WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8.root"},
	                               {"ST_tW_top_5f", "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"ST_tW_antitop_5f", "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"ST_t-channel_top_4f", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8.root"},
	                               {"ST_t-channel_antitop_4f", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8.root"},
	                           	   {"ST_t-channel_antitop_5f", "ST_t-channel_antitop_5f_TuneCP5_13TeV-powheg-pythia8.root"},
	                               {"ST_t-channel_top_5f", "ST_t-channel_top_5f_TuneCP5_13TeV-powheg-pythia8.root"}};

	map<TString, float> sub16XSEC = {{"DY",1208 },
	                               {"WJetsToQQ_HT180",3105 },
	                               {"ST_tW_top_5f", 38.09},
	                               {"ST_tW_antitop_5f",38.09},
	                               {"ST_t-channel_top_4f", 35.6},
	                               {"ST_t-channel_antitop_4f",35.6}};

    map<TString, float> sub17XSEC = {{"DY", 1728},
	                               {"WJetsToQQ_HT400",1447 },
	                               {"WJetsToQQ_HT600",318.8 },
	                               {"ST_tW_top_5f", 34.91},
	                               {"ST_tW_antitop_5f",34.97},
	                               {"ST_t-channel_top_4f", 113.3},
	                               {"ST_t-channel_antitop_4f",67.91}};

    map<TString, float> sub18XSEC = {{"DY", 1728},
	                               {"WJetsToQQ_HT400",1447 },
	                               {"WJetsToQQ_HT600",318.8 },
	                               {"ST_tW_top_5f", 34.91},
	                               {"ST_tW_antitop_5f",34.97},
	                               {"ST_t-channel_top_4f", 113.3},
	                               {"ST_t-channel_antitop_4f",67.91},
	                           	   {"ST_t-channel_antitop_5f",71.74},
	                           	   {"ST_t-channel_top_5f",119.7}};

	subdominantBkgFiles.insert(pair<TString, map<TString, TString>>("2016",sub16));
	subdominantBkgFiles.insert(pair<TString, map<TString, TString>>("2017",sub17));
	subdominantBkgFiles.insert(pair<TString, map<TString, TString>>("2018",sub18));

	subdominantBkgXSEC.insert(pair<TString, map<TString, float>>("2016",sub16XSEC));
	subdominantBkgXSEC.insert(pair<TString, map<TString, float>>("2017",sub17XSEC));
	subdominantBkgXSEC.insert(pair<TString, map<TString, float>>("2018",sub18XSEC));
//----------------------------------------------------------------------------------------------------------------	
	//Nominal TT files (MC):
	map<TString, TString> eosNomTT16 = {{"TTNominal", "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"}};
	
    map<TString, TString> eosNomTT17 = {{"TTHadronic_0", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	   {"TTSemiLeptonic_0", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	   {"TTTo2L2Nu_0", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};	                              

    map<TString, TString> eosNomTT18 = {{"TTHadronic_0", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	   {"TTSemiLeptonic_0", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"},
	                               	   {"TTTo2L2Nu_0", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"}};
	                                  
	map<TString, float> eosNomTTXSEC16 = {{"TTNominal",832.}};

    map<TString, float> eosNomTTXSEC17 = {{"TTHadronic_0",377.96},
	                                    {"TTSemiLeptonic_0",365.34},
	                                    {"TTTo2L2Nu_0",88.29}};	                              

    map<TString, float> eosNomTTXSEC18 = {{"TTHadronic_0",377.96},
	                                    {"TTSemiLeptonic_0",365.34},
	                                    {"TTTo2L2Nu_0",88.29}};		                                  

    ttNominalFiles.insert(pair<TString, map<TString, TString>>("2016", eosNomTT16));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2017", eosNomTT17));
	ttNominalFiles.insert(pair<TString, map<TString, TString>>("2018", eosNomTT18));                               

	ttNominalXSEC.insert(pair<TString, map<TString, float>>("2016",eosNomTTXSEC16));
	ttNominalXSEC.insert(pair<TString, map<TString, float>>("2017",eosNomTTXSEC17));
	ttNominalXSEC.insert(pair<TString, map<TString, float>>("2018",eosNomTTXSEC18));

//----------------------------------------------------------------------------------------------------------------


	//thsese will be used either for naming the files for FillHistograms, or either for getting the histo names 
	map<TString, TString> files2016 = {{"data",  "2016/Histo_Data_2016_100.root"},
	                                   {"mcSig", "2016/Histo_TT_Mtt-700toInf_100.root"},
	                                   {"mcSub", "2016/Histo_SubdominantBkgs_100.root"},
	                               	   {"qcd"  , "2016/Histo_QCD_HT300toInf_100.root"}};

	map<TString, TString> files2016Loose = {{"data",  "2016/Histo_Data_2016_100_Loose.root"},
	                                   		{"mcSig", "2016/Histo_TT_Mtt-700toInf_100_Loose.root"},
	                                   		{"mcSub", "2016/Histo_SubdominantBkgs_100_Loose.root"},
	                                   		{"qcd"  , "2016/Histo_QCD_HT300toInf_100_Loose.root"}};	                                   
																		 
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
	                                   																			 
	files.insert(pair<TString, map<TString, TString>>("2016", files2016));
	files.insert(pair<TString, map<TString, TString>>("2017", files2017));
	files.insert(pair<TString, map<TString, TString>>("2018", files2018));
	files.insert(pair<TString, map<TString, TString>>("2016_Loose", files2016Loose));
	files.insert(pair<TString, map<TString, TString>>("2017_Loose", files2017Loose));
	files.insert(pair<TString, map<TString, TString>>("2018_Loose", files2018Loose));


	map<TString, TString> files2016Reduced = {{"data",  "2016/Histo_Data_2016_100_reduced.root"},
	                                   		  {"mcSig", "2016/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2016/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2016/Histo_QCD_HT300toInf_100_reduced.root"}};
																		 
	map<TString, TString> files2017Reduced = {{"data",  "2017/Histo_Data_2017_100_reduced.root"},
	                                   		  {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "20167Histo_QCD_HT300toInf_100_reduced.root"}};                                		  
	
	map<TString, TString> files2018Reduced = {{"data",  "2018/Histo_Data_2018_100_reduced.root"},
	                                   	      {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced.root"},
	                                   		  {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced.root"},
	                                   		  {"qcd"  , "2018/Histo_QCD_HT300toInf_100_reduced.root"}};

	map<TString, TString> files2016ReducedLoose = {{"data",  "2016/Histo_Data_2016_100_reduced_Loose.root"},
	                                   		       {"mcSig", "2016/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		       {"mcSub", "2016/Histo_SubdominantBkgs_100_reduced_Loose.root"},
	                                   		       {"qcd"  , "2016/Histo_QCD_HT300toInf_100_reduced_Loose.root"}};
																		 
	map<TString, TString> files2017ReducedLoose = {{"data",  "2017/Histo_Data_2017_100_reduced_Loose.root"},
	                                   		  	   {"mcSig", "2017/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		  	   {"mcSub", "2017/Histo_SubdominantBkgs_100_reduced_Loose.root"},
	                                   		  	   {"qcd"  , "2017/Histo_QCD_HT300toInf_100_reduced_Loose.root"}};                                		  
	
	map<TString, TString> files2018ReducedLoose = {{"data",  "2018/Histo_Data_2018_100_reduced_Loose.root"},
	                                   	           {"mcSig", "2018/Histo_TT_Mtt-700toInf_100_reduced_Loose.root"},
	                                   		       {"mcSub", "2018/Histo_SubdominantBkgs_100_reduced_Loose.root"},
	                                   		       {"qcd"  , "2018/Histo_QCD_HT300toInf_100_reduced_Loose.root"}};                                  		  
																		 
	filesReduced.insert(pair<TString, map<TString, TString>>("2016", files2016Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017", files2017Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018", files2018Reduced));
	filesReduced.insert(pair<TString, map<TString, TString>>("2016Loose", files2016ReducedLoose));
	filesReduced.insert(pair<TString, map<TString, TString>>("2017Loose", files2017ReducedLoose));
	filesReduced.insert(pair<TString, map<TString, TString>>("2018Loose", files2018ReducedLoose));
	
	
	luminosity.insert(pair<TString, float>("2016",35920));
	luminosity.insert(pair<TString, float>("2017",41530));
	luminosity.insert(pair<TString, float>("2018",59740));

	luminosityCR.insert(pair<TString, float>("2016",1670));
	luminosityCR.insert(pair<TString, float>("2017",41530));
	luminosityCR.insert(pair<TString, float>("2018",59740));

	if(!isLoose)
	{
		deepCSVFloatMap.insert(pair<TString, float>("2016",0.6321));
		deepCSVFloatMap.insert(pair<TString, float>("2017",0.4941));
		deepCSVFloatMap.insert(pair<TString, float>("2018",0.4184));

		topTaggerCuts.insert(pair<TString, float>("2016",0.2));
		topTaggerCuts.insert(pair<TString, float>("2017",0.0));
		topTaggerCuts.insert(pair<TString, float>("2018",0.1));
	}
	else
	{
		deepCSVFloatMap.insert(pair<TString, float>("2016",0.2217));
		deepCSVFloatMap.insert(pair<TString, float>("2017",0.1522));
		deepCSVFloatMap.insert(pair<TString, float>("2018",0.1241));

		topTaggerCuts.insert(pair<TString, float>("2016",-0.3));
		topTaggerCuts.insert(pair<TString, float>("2017",-0.4));
		topTaggerCuts.insert(pair<TString, float>("2018",-0.4));
	}

	triggerSRConst.insert(pair<TString, int>("2016",2));
	triggerSRConst.insert(pair<TString, int>("2017",5));
	triggerSRConst.insert(pair<TString, int>("2018",5));

	triggerCRConst.insert(pair<TString, int>("2016",4));
	triggerCRConst.insert(pair<TString, int>("2017",5));
	triggerCRConst.insert(pair<TString, int>("2018",5));

	variableConstant.insert(pair<TString, int>("mJJ",  0));
	variableConstant.insert(pair<TString, int>("ptJJ", 1));
	variableConstant.insert(pair<TString, int>("yJJ",  2));
	variableConstant.insert(pair<TString, int>("jetPt0", 3));
	variableConstant.insert(pair<TString, int>("jetPt1", 4));
	variableConstant.insert(pair<TString, int>("jetY0", 5));
	variableConstant.insert(pair<TString, int>("jetY1", 6));

	//these are fit results from mixed medium wp and loose wp
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2016", 3.2755e+03));
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2017", 3.3491e+03));
	Nbkg2Constants.insert(pair<TString, float>("Nbkg2018", 4.5360e+03));

	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2016_error", 2.09e+02));
	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2017_error", 6.98e+02));
	Nbkg2ConstantsErrors.insert(pair<TString, float>("Nbkg2018_error", 1.52e+02));
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
