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
             "2018":['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_20UL.root',
                     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_20UL.root',
                     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_20UL.root']}

eospath = '/eos/cms/store/user/gbakas/ttbar/topTagger/ul-{}/Signal/'.format(year)


for ifile in  all_files[year]:
    split_file_underscore = ifile.split('_')
    print(split_file_underscore[0])
    print(ifile)
    for mjj_cut in mJJCuts:
        print('Current mjj cut:', mjj_cut)
        #os.system(f'root -l -b -q \'FillHistograms_Extended_PS_PDF.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{weightType}\",\"{mjj_cut}\")\'')
        os.system(f'root -l -b -q \'FillHistograms_Reduced_PS_PDF.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{weightType}\",\"{mjj_cut}\")\'')
