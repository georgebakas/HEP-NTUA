import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
weightType = sys.argv[2]
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


for ifile in  all_files[year]:
    split_file_underscore = ifile.split('_')
    print(split_file_underscore[0])
    print(ifile)
    os.system(f'root -l -b -q \'ResponseMatrices_PS_PDF.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{weightType}\")\'')
    #os.system(f'root -l -b -q \'FillHistograms_Reduced_PS_PDF.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{weightType}\")\'')
