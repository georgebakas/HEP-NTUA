import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')
print(f'Type of file you wan to add: {str(sys.argv[2])}')

year = sys.argv[1]
weightType = sys.argv[2]

tof = 'Reduced'

allFiles = []

ps_weights = {2:"isrRedHi", 3:"fsrRedHi", 4:"isrRedLo", 5:"fsrRedLo", 6:"isrDefHi", 7:"fsrDefHi",
              8:"isrDefLo", 9:"fsrDefLo", 10:"isrConHi", 11:"fsrConHi", 12:"isrConLo", 13:"fsrConLo"}

for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/HistoReduced_TTToHadronic*.root'.format(year,weightType,tof), recursive=True)):
#for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/*.root'.format(year,weightType,tof), recursive=True)): #for responses
    #print('file: {}'.format(file_name))
    split_file_name = file_name.split('/')
    split_file_underscore = split_file_name[-1].split('_')

    print(split_file_underscore[1])
    tt_process_name = split_file_underscore[1]
    if weightType != 'PSWeights':
        dot_split = split_file_underscore[-1].split('.')
        weight_sufix = dot_split[0]
        print(weight_sufix)
        #print(split_file_underscore[2])
        if weightType == 'PDFWeights':
             pdf_or_scale = 'pdf'
        else:
            pdf_or_scale = 'scale'
        #print(f'hadd -f {year}/{weightType}/combined/Histo{tof}_TT_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_TTToHadronic_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_TTToSemiLeptonic_{pdf_or_scale}_{weight_sufix}.root {year}/{weightType}/Histo{tof}_TTTo2L2Nu_{pdf_or_scale}_{weight_sufix}.root')
        os.system(f'hadd -f {year}/Responses{weightType}/combined/Histo{tof}_TT_{pdf_or_scale}_{weight_sufix}.root {year}/Responses{weightType}/Histo{tof}_TTToHadronic_{pdf_or_scale}_{weight_sufix}.root {year}/Responses{weightType}/Histo{tof}_TTToSemiLeptonic_{pdf_or_scale}_{weight_sufix}.root {year}/Responses{weightType}/Histo{tof}_TTTo2L2Nu_{pdf_or_scale}_{weight_sufix}.root')
    else:

        dot_split = split_file_underscore[2].split('.')
        weight_sufix = dot_split[0]
        print(weight_sufix)
        os.system(f'hadd -f {year}/Responses{weightType}/combined/Histo{tof}_TT_{weightType}_{weight_sufix}.root {year}/Responses{weightType}/Histo{tof}_TTToHadronic_{weight_sufix}.root {year}/Responses{weightType}/Histo{tof}_TTToSemiLeptonic_{weight_sufix}.root {year}/Responses{weightType}/Histo{tof}_TTTo2L2Nu_{weight_sufix}.root')
