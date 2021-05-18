import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
weightType = sys.argv[2]
allFiles = []

all_files = {"2016_preVFP":['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_20UL.root',
                     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_20UL.root',
                     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_20UL.root '],
             "2016_postVFP":['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_20UL.root',
                     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_20UL.root',
                     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_20UL.root'],
             "2017":['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_19UL.root',
                     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_19UL.root',
                     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_19UL.root'],
             "2018":['TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root',
                     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root',
                     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root']}

if year!='2018':
        eospath = '/eos/cms/store/user/gbakas/ttbar/topTagger/ul-{}/Signal/'.format(year)
else:
        eospath = '/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2018/'


for ifile in  all_files[year]:
    split_file_underscore = ifile.split('_')
    print(split_file_underscore[0])
    print(ifile)
    os.system(f'root -l -b -q \'ResponseMatrices_PS_PDF.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{weightType}\")\'')
