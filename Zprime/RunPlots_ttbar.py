import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]

all_files = {"2016":['TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root',
                     'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root',
                     'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root'],
             "2017":['TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root',
                     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root',
                     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root'],
             "2018":['TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root',
                     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root',
                     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root']}

eospath = '/eos/cms/store/user/gbakas/ttbar/topTagger/mc-{}/Signal/'.format(year)


for ifile in  all_files[year]:
    split_file_underscore = ifile.split('_')
    print(split_file_underscore[0])
    print(ifile)
    os.system(f'root -l -q \'PlotCutFlow_JetMassSoftDrop.C(\"{year}\", \"HistoCutFlowJetMassSoftDrop_{ifile}\",\"{split_file_underscore[0]}\")\'')
    #break
