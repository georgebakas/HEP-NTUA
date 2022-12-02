import os 

years = ["2016_preVFP", "2016_postVFP", "2017", "2018"]


for year in  years:
    for i in range(4):
        if i == 1: 
            continue
        os.system(f'root -l -b -q \'TagAndProbe.C(\"{year}\",{i})\'')