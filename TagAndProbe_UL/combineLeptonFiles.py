import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
weightType = 'Nominal'

for ifile, file_name in enumerate(glob.iglob('{}/{}/LeptonAnalysis_TTToHadronic*.root'.format(year,weightType), recursive=True)):
    print('file: {}'.format(file_name))
    split_file_name = file_name.split('/')
    split_file_underscore = split_file_name[-1].split('_')

    tt_process_name = split_file_underscore[1]
    print(ifile, file_name)
    os.system(f'hadd -f {year}/{weightType}/combined/LeptonAnalysis_Nominal.root {year}/{weightType}/LeptonAnalysis_TTToHadronic.root {year}/{weightType}/LeptonAnalysis_TTToSemiLeptonic.root {year}/{weightType}/LeptonAnalysis_TTTo2L2Nu.root')

