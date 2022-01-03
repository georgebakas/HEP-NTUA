import os
import sys

years = ["2016_preVFP", "2016_postVFP", "2017", "2018"]
btagvar = ["CR", "SR"]
for year in years:
    os.system(f'rm {year}/StackPlots/*.pdf')
    for region in btagvar:
        print('year', year)
        os.system(f'root -l -b -q \'plotStackHisto.C(\"{year}\",1000, \"{region}\")\'')
        #break
    #break