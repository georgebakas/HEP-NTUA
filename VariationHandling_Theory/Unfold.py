import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')
print(f'Type of file you wan to add: {str(sys.argv[2])}')

year = sys.argv[1]
weightType = sys.argv[2]

allFiles = []

ps_weights = {2:"isrRedHi", 3:"fsrRedHi", 4:"isrRedLo", 5:"fsrRedLo", 6:"isrDefHi", 7:"fsrDefHi",
              8:"isrDefLo", 9:"fsrDefLo", 10:"isrConHi", 11:"fsrConHi", 12:"isrConLo", 13:"fsrConLo"}

for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/*TTToHadronic*.root'.format(year,weightType), recursive=True)):
    print(file_name)
    if weightType != 'SystematicsFiles':
        split_file_name = file_name.split('/')
        split_file_underscore = split_file_name[-1].split('_')

        print(split_file_underscore[1])
        tt_process_name = split_file_underscore[1]

        dot_split = split_file_underscore[-1].split('.')
        weight_sufix = dot_split[0]
        print(weight_sufix)

        #os.system(f'root -l -b -q \'Unfold_data.cpp(\"{year}\",\"{weightType}\",\"{weight_sufix}\", true)\'')
    elif weightType == 'Nominal':
        weight_sufix = ""
    else:
        #TT, and TTJets files are handled below
        split_file_name = file_name.split('/')
        split_file_underscore = split_file_name[-1].split('_')

        print(split_file_underscore[1])
        tt_process_name = split_file_underscore[1]

        dot_split = split_file_name[-1].replace(split_file_underscore[0]+'_TTToHadronic_', '')
        dot_split = dot_split.split('.')
        weight_sufix = dot_split[0]
        print(weight_sufix)
        #(TString inYear, TString dir, TString inputFile, bool isThreeProcesses, bool isParton = true, int unfoldMethod=1)

    os.system(f'root -l -b -q \'Unfold_data.cpp(\"{year}\",\"{weightType}\",\"{weight_sufix}\", true)\'')
    #break

combined_files = ['TT','TTJets']
#check for files of type TT_Tune
if weightType == 'SystematicsFiles':
    for tt in combined_files:
        for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/*{}_*.root'.format(year,weightType, tt), recursive=True)):
            split_file_name = file_name.split('/')
            split_file_underscore = split_file_name[-1].split('_')

            print(split_file_underscore[1])
            tt_process_name = split_file_underscore[1]

            dot_split = split_file_name[-1].replace('ResponsesEfficiency_', '')
            dot_split = dot_split.split('.')
            weight_sufix = dot_split[0]
            print(weight_sufix)
            os.system(f'root -l -b -q \'Unfold_data.cpp(\"{year}\",\"{weightType}\",\"{weight_sufix}\", false)\'')
#eof
