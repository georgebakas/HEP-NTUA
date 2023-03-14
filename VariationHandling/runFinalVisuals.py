import sys
import os 

years = ["2016_preVFP", "2016_postVFP", "2017", "2018"]

for year in years:
    os.system(f'root -l -b -q \'PlotMassDiff.C(\"{year}\")\'')

# Stack histograms
# void plotStackHisto(TString year, int mJJCut = 1000, TString region = "SR")

# Final Results Fiducial
norms = ['false', 'true']
for norm in norms:
    os.system(f'root -l -b -q \'FinalResults_Fiducial.cxx(\"{norm}\")\'')
    os.system(f'root -l -b -q \'Fiducial_Basic.C(\"{norm}\")\'')

phase_space = ['Parton', 'Particle']

for p in phase_space:

    # Final Results 
    norm = 'false'
    os.system(f'root -l -b -q \'FinalResults.cxx(\"{p}\",{norm})\'')
    norm = 'true'
    os.system(f'root -l -b -q \'FinalResults.cxx(\"{p}\",{norm})\'')

    # Chi Square evaluation and results 
    os.system(f'root -l -b -q \'ExportChi2.C(\"{p}\")\'')
    os.system(f'root -l -b -q \'ExportChi2Norm.C(\"{p}\")\'')

    # Final Results 
    os.system(f'root -l -b -q \'DrawWithChiSquare.cxx(\"{p}\")\'')
    os.system(f'root -l -b -q \'DrawWithChiSquare_normalised.cxx(\"{p}\")\'')