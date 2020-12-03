import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

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

variations = ['up', 'down']


for var in variations:
    for ifile in  all_files[year]:
        print('variation is:', var)
        split_file_underscore = ifile.split('_')
        print(split_file_underscore[0])
        print(ifile)
        os.system(f'root -l -b -q \'FillHistograms_Extended_bTag.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{var}\")\'')
        #os.system(f'root -l -b -q \'FillHistograms_Reduced_bTag.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{var}\")\'')
       #break
    #break
