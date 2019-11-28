#signal
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8", true)	
#mtt files
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8",true)
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8",true)

#QCD bkg
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true)
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true)
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true)
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true)
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true)
root -b -l -q CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true)