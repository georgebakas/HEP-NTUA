import sys
import os 

years = ["2016_preVFP", "2016_postVFP", "2017", "2018"]

phase_space = ['Parton', 'Particle']

for y in years:
    for region in ['SR', 'CR']:

        # Plot Stack histograms without SF without Magenta
        os.system(f'root -l -b -q \'plotStackHisto.C(\"{y}\", 1000, \"{region}\", false)\'')

        # Plot Stack histograms with SF without Magenta
        os.system(f'root -l -b -q \'plotStackHisto.C(\"{y}\", 1000, \"{region}\", true)\'')

        # Plot Stack histograms with SF with Magenta
        os.system(f'root -l -b -q \'plotStackHisto_TopTaggerErrors.C(\"{y}\", 1000, \"{region}\", true)\'')