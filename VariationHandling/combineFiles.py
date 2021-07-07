import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')
print(f'Type of file you want to add: {str(sys.argv[2])}')

year = sys.argv[1]
weightType = sys.argv[2]
reduced = sys.argv[3] #type of file. Histo or HistoReduced
print(reduced)
if reduced == 'True':
    tof = 'Reduced'
else:
    tof = ''
allFiles = []

ps_weights = {2:"isrRedHi", 3:"fsrRedHi", 4:"isrRedLo", 5:"fsrRedLo", 6:"isrDefHi", 7:"fsrDefHi",
              8:"isrDefLo", 9:"fsrDefLo", 10:"isrConHi", 11:"fsrConHi", 12:"isrConLo", 13:"fsrConLo"}

mJJCuts = [1000, 1200, 1400, 1600, 1800, 2000]

for mJJCut in mJJCuts:
    for ifile, file_name in enumerate(glob.iglob('{}/{}/Histo{}_{}_TTToHadronic*.root'.format(year,weightType,tof,mJJCut), recursive=True)):
    #for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/TTToHadronic*.root'.format(year,weightType,tof), recursive=True)): #for responses
        print('file: {}'.format(file_name))
        split_file_name = file_name.split('/')
        split_file_underscore = split_file_name[-1].split('_')

        tt_process_name = split_file_underscore[1]
        if weightType not in  ['PSWeights', 'bTagVariation', 'Nominal', 'JES']:
            dot_split = split_file_underscore[-1].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            #print(split_file_underscore[2])
            if weightType == 'PDFWeights':
                 pdf_or_scale = 'pdf'
            else:
                pdf_or_scale = 'scale'
            #print(f'hadd -f {year}/{weightType}/combined/Histo{tof}_{mJJCut}_TT_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_TTToHadronic_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_TTToSemiLeptonic_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_TTTo2L2Nu_{pdf_or_scale}_{weight_sufix}.root')
            os.system(f'hadd -f {year}/{weightType}/combined/Histo{tof}_{mJJCut}_TT_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToHadronic_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToSemiLeptonic_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTTo2L2Nu_{pdf_or_scale}_{weight_sufix}.root')
        elif weightType == 'PSWeights':
            dot_split = split_file_underscore[3].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'hadd -f {year}/{weightType}/combined/Histo{tof}_{mJJCut}_TT_{weightType}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToHadronic_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToSemiLeptonic_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTTo2L2Nu_{weight_sufix}.root')

        elif weightType == 'bTagVariation':
            print(ifile, file_name)
            dot_split = split_file_underscore[3].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'hadd -f {year}/{weightType}/combined/Histo{tof}_{mJJCut}_TT_{weightType}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToHadronic_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToSemiLeptonic_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTTo2L2Nu_{weight_sufix}.root')

        elif weightType == 'JES':
            print(ifile, file_name)
            dot_split = split_file_underscore[3].split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'hadd -f {year}/{weightType}/combined/Histo{tof}_{mJJCut}_TT_{weightType}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToHadronic_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToSemiLeptonic_{weight_sufix}.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTTo2L2Nu_{weight_sufix}.root')

        else:
            print(ifile, file_name)
            os.system(f'hadd -f {year}/{weightType}/combined/Histo{tof}_{mJJCut}_TT.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToHadronic.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTToSemiLeptonic.root {year}/{weightType}/Histo{tof}_{mJJCut}_TTTo2L2Nu.root')
    break
