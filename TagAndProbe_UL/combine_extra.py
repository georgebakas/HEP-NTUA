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
    for ifile, file_name in enumerate(glob.iglob('{}/{}/TTbarExtraAnalysis_TTToHadronic*.root'.format(year,weightType), recursive=True)):
    #for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/TTToHadronic*.root'.format(year,weightType,tof), recursive=True)): #for responses
        print('file: {}'.format(file_name))
        split_file_name = file_name.split('/')
        split_file_underscore = split_file_name[-1].split('_')
        tt_process_name = split_file_underscore[1]
        print(ifile, file_name)
        os.system(f'hadd -f {year}/{weightType}/combined/TTbarExtraAnalysis_TT_Nominal.root {year}/{weightType}/TTbarExtraAnalysis_TTToHadronic.root {year}/{weightType}/TTbarExtraAnalysis_TTToSemiLeptonic.root {year}/{weightType}/TTbarExtraAnalysis_TTTo2L2Nu.root')
        
        
    break
