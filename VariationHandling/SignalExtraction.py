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

for ifile, file_name in enumerate(glob.iglob('{}/{}/MassFitResults_*.root'.format(year,weightType), recursive=True)):

    print(ifile, file_name)
    split_file_name = file_name.split('/')  
    os.system(f'root -l -b -q \'SignalExtraction_UnequalBins.cpp(\"{year}\",\"{weightType}\",\"{split_file_name[-1]}\")\'')
