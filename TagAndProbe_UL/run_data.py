import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

all_files = {"2016_preVFP":['JetHT_Run2016-21Feb2020_UL2016_HIPM-v1.root'],
             "2016_postVFP":['JetHT_Run2016-21Feb2020_UL2016-v1.root'],
             "2017":['JetHT_Run2017-UL2017_MiniAODv2-v1.root'],
             "2018":['JetHT_Run2018-UL2018_MiniAODv2-v1.root']}

eospath = '/eos/cms/store/user/gbakas/ttbar/JetHT/ul-{}/'.format(year)
mJJCuts = [1000]#, 1200, 1400, 1600, 1800, 2000]

for ifile in  all_files[year]:
    split_file_underscore = ifile.split('_')
    print(split_file_underscore[0])
    print(ifile)
    for mjj_cut in mJJCuts:
        print('Current mjj cut:', mjj_cut)
        os.system(f'root -l -b -q \'TagAndProbe_extra_data.C(\"{eospath+ifile}\",\"{year}\",{mjj_cut})\'')