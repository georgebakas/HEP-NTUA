import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

all_files = {"2016_preVFP":['TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root'],
             "2016_postVFP":['TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root'],
             "2017":['TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root'],
             "2018":['TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root']}

eospath = '/eos/cms/store/user/ipapakri/ttbar/MC/Signal/{}/'.format(year)

variations = ['up', 'down']

for var in variations:
    for ifile in  all_files[year]:
        print('variation is:', var)
        split_file_underscore = ifile.split('_')
        print(split_file_underscore[0])
        print(ifile)
        os.system(f'root -l -b -q \'GetHistograms_bTag.C(\"{eospath+ifile}\",\"{split_file_underscore[0]}\",\"{year}\",\"{var}\")\'')
        #break
    #break
