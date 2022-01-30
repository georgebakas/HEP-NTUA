import sys
import os

# I need to map these with the correct mass Cuts
mass_cut = {"1200_12":1000, "1400_14":1200, "1600_16":1400, "1800_18":1600, "2000_20":1600,
            "2500_25":2000, "3000_30":2000, "3500_35":2000, "4000_40":2000, "4500_45":2000}


for key, value in mass_cut.items():
    datacard = f'datacard_chi_SR_mZ_{key}_cut_{value}.txt'
    if key != "1400_14":
        cmd = f'combineCards.py ../2016_preVFP_test/{datacard} ../2017_test/{datacard} ../2018_test/{datacard} > {datacard}'
    else:
        cmd = f'combineCards.py ../2017_test/{datacard} ../2018_test/{datacard} > {datacard}'
    os.system(cmd)