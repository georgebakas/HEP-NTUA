import os 
import sys

years = ["2016_preVFP", "2016_postVFP", "2017", "2018"]

# for year in  years:

year = sys.argv[1]
for i in range(1,4):
    if i != 1: 
        os.system(f'root -l -b -q \'TopTaggerDatavsMCOutput.C(\"{year}\",{i})\'')
        os.system(f'root -l -b -q \'TopTaggerDataVSMCOutput_pTRegions.C(\"{year}\",{i})\'')