import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')
print(f'Type of file you want to add: {str(sys.argv[2])}')

year = sys.argv[1]
weightType = sys.argv[2]

allFiles = []

ps_weights = {2:"isrRedHi", 3:"fsrRedHi", 4:"isrRedLo", 5:"fsrRedLo", 6:"isrDefHi", 7:"fsrDefHi",
            8:"isrDefLo", 9:"fsrDefLo", 10:"isrConHi", 11:"fsrConHi", 12:"isrConLo", 13:"fsrConLo"}

mJJCuts = [1000] #, 1200, 1400, 1600, 1800, 2000]
vars_ = ['PSWeights', 'bTagVariation', 'Nominal', 'JES', 'SystematicsFiles']

for mJJCut in mJJCuts:
    for ifile, file_name in enumerate(glob.iglob('{}/{}/TagAndProbeHisto_TT_{}_TTToHadronic*_newBins.root'.format(year,weightType,mJJCut), recursive=True)):
    #for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/TTToHadronic*_newBins.root'.format(year,weightType,tof), recursive=True)): #for responses
        print('file: {}'.format(file_name))
        split_file_name = file_name.split('/')
        split_file_underscore = split_file_name[-1].split('_')

        tt_process_name = split_file_underscore[1]
        if weightType not in vars_:
            dot_split = split_file_underscore[-1].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            #print(split_file_underscore[2])
            if weightType == 'PDFWeights':
                pdf_or_scale = 'pdf'
            else:
                pdf_or_scale = 'scale'
            os.system(f'hadd -f {year}/{weightType}/combined/TagAndProbeHisto_{mJJCut}_TT_{pdf_or_scale}_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToHadronic_{pdf_or_scale}_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToSemiLeptonic_{pdf_or_scale}_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTTo2L2Nu_{pdf_or_scale}_{weight_sufix}_newBins.root')
        
        elif weightType == 'PSWeights':
            dot_split = split_file_underscore[3].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'hadd -f {year}/{weightType}/combined/TagAndProbeHisto_{mJJCut}_TT_{weightType}_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToHadronic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToSemiLeptonic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTTo2L2Nu_{weight_sufix}_newBins.root')

        elif weightType == 'bTagVariation':
            print(ifile, file_name)
            dot_split = split_file_underscore[3].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'hadd -f {year}/{weightType}/combined/TagAndProbeHisto_{mJJCut}_TT_{weightType}_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToHadronic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToSemiLeptonic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTTo2L2Nu_{weight_sufix}_newBins.root')

        elif weightType == 'JES':
            print(ifile, file_name)
            dot_split = split_file_underscore[3].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'hadd -f {year}/{weightType}/combined/TagAndProbeHisto_{mJJCut}_TT_{weightType}_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToHadronic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToSemiLeptonic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTTo2L2Nu_{weight_sufix}_newBins.root')
        
        elif weightType == 'SystematicsFiles':
            print(ifile, file_name)
            dot_split = split_file_underscore[3].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'hadd -f {year}/{weightType}/combined/TagAndProbeHisto_{mJJCut}_TT_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToHadronic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToSemiLeptonic_{weight_sufix}_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTTo2L2Nu_{weight_sufix}_newBins.root')
        
        else:
            print(ifile, file_name)
            os.system(f'hadd -f {year}/{weightType}/combined/TagAndProbeHisto_{mJJCut}_TT_Nominal_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToHadronic_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTToSemiLeptonic_newBins.root {year}/{weightType}/TagAndProbeHisto_TT_{mJJCut}_TTTo2L2Nu_newBins.root')
        
        
    break
