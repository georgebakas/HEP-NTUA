import os 

years = ["2016_preVFP", "2016_postVFP", "2017", "2018"]


for year in  years:
    for i in range(4):
        if i == 1: 
            continue
        if i == 3:
            #os.system(f'root -l -b -q \'FillHistograms.C(\"{year}\",{i})\'')
            os.system(f'root -l -b -q \'FillHistograms_Reduced_UnequalBinning.cpp(\"{year}\",{i})\'')