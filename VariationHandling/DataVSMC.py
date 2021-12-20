import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
isExtended = sys.argv[2]
allFiles = []

all_files = {"2016_preVFP":['TTToHadronic',
                     'TTToSemiLeptonic',
                     'TTTo2L2Nu',
                     'ST_tW_top_5f_NoFullyHadronicDecays',
                     'ST_tW_antitop_5f_NoFullyHadronicDecays',
                     'ST_t-channel_top_5f_InclusiveDecays',
                     'ST_t-channel_top_4f_InclusiveDecays',
                     'ST_t-channel_antitop_5f_InclusiveDecays',
                     'QCD_HT700to1000',
                     'QCD_HT500to700',
                     'QCD_HT300to500',
                     'QCD_HT2000toInf',
                     'QCD_HT1500to2000',
                     'QCD_HT1000to1500'],

             "2016_postVFP":['TTToHadronic',
                     'TTToSemiLeptonic',
                     'TTTo2L2Nu',
                     'ST_tW_top_5f_NoFullyHadronicDecays',
                     'ST_tW_antitop_5f_NoFullyHadronicDecays',
                     'ST_t-channel_top_5f_InclusiveDecays',
                     'ST_t-channel_top_4f_InclusiveDecays',
                     'ST_t-channel_antitop_5f_InclusiveDecays',
                     'QCD_HT700to1000',
                     'QCD_HT500to700',
                     'QCD_HT300to500',
                     'QCD_HT2000toInf',
                     'QCD_HT1500to2000',
                     'QCD_HT1000to1500'],

             "2017":['TTToHadronic',
                     'TTToSemiLeptonic',
                     'TTTo2L2Nu',
                     'QCD_HT1000to1500',
                     'ST_tW_top_5f_NoFullyHadronicDecays',
                     'ST_tW_antitop_5f_NoFullyHadronicDecays',
                     'ST_t-channel_top_5f_InclusiveDecays',
                     'ST_t-channel_top_4f_InclusiveDecays',
                     'ST_t-channel_antitop_5f_InclusiveDecays',
                     'QCD_HT700to1000',
                     'QCD_HT500to700',
                     'QCD_HT300to500',
                     'QCD_HT2000toInf',
                     'QCD_HT1500to2000'],
             
             "2018":['TTToHadronic',
                     'TTToSemiLeptonic',
                     'TTTo2L2Nu',
                     'ST_tW_top_5f_NoFullyHadronicDecays',
                     'ST_tW_antitop_5f_NoFullyHadronicDecays',
                     'ST_t-channel_top_5f_InclusiveDecays',
                     'ST_t-channel_top_4f_InclusiveDecays',
                     'ST_t-channel_antitop_5f_InclusiveDecays',
                     'ST_t-channel_antitop_4f_InclusiveDecays',
                     'QCD_HT700to1000',
                     'QCD_HT500to700',
                     'QCD_HT300to500',
                     'QCD_HT2000toInf',
                     'QCD_HT1500to2000',
                     'QCD_HT1000to1500']}

eospath = ['/eos/cms/store/user/gbakas/ttbar/topTagger/ul-{}/'.format(year),
            '/eos/cms/store/user/gbakas/ttbar/JetHT/ul-{}/'.format(year)]


all_files_data = {
    "2016_preVFP":"JetHT_Run2016-21Feb2020_UL2016_HIPM-v1.root",
    "2016_postVFP":"JetHT_Run2016-21Feb2020_UL2016-v1.root",
    "2017":"JetHT_Run2017-UL2017_MiniAODv2-v1.root",
    "2018":"JetHT_Run2018-UL2018_MiniAODv2-v1.root"
}

mJJCuts = [1000]#, 1200, 1400, 1600, 1800, 2000]

for ifile in  all_files[year]:
    print(ifile)
    if isExtended == 'True':
            os.system(f'root -l -b -q \'DataVSMC_Extended.C(\"{eospath[0]}", \"{ifile}\", \"false\", \"{year}\")\'')
    else:
            os.system(f'root -l -b -q \'DataVSMC_Reduced.C(\"{eospath[0]}", \"{ifile}\", \"false\", \"{year}\")\'')

#if isExtended == 'True':
#    split_file_underscore = all_files_data[year].split('_')
#    os.system(f'root -l -b -q \'DataVSMC_Extended.C(\"{eospath[1]}", \"{all_files_data[year]}\", \"true\", \"{year}\")\'')
#else:
#    os.system(f'root -l -b -q \'DataVSMC_Extended.C(\"{eospath[1]}", \"{all_files_data[year]}\", \"true\", \"{year}\")\'')