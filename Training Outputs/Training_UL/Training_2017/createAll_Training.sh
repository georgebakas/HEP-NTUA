#signal
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Signal/TTToHadronic_TuneCP5_13TeV-powheg-pythia8_19UL", true)'
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Signal/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_19UL", true)'
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Signal/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_19UL", true)'	

#QCD bkg
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Bkg/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL",false)'
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Bkg/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL",false)'
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Bkg/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL",false)'
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Bkg/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL",false)'
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Bkg/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL",false)'
root -b -l -q 'CreateTreeForMVATraining_Matching.C("/eos/cms/store/user/gbakas/ttbar/topTagger/ul-2017/Bkg/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8_19UL",false)'