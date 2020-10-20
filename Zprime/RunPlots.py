import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
level = sys.argv[2]
mJJCut = sys.argv[3]
allFiles = []
histo_names = []

for ifile, file_name in enumerate(glob.iglob('{}_backup/HistoMassWindows*.root'.format(year), recursive=True)):
    print('file: {}'.format(file_name))
    file_to_append = file_name.split('/')
    allFiles.append(file_to_append[-1])
    split_file_in_underscores = file_to_append[-1].split('_')
    mass_name = split_file_in_underscores[2]
    width_name = split_file_in_underscores[3]
    histo_names.append(mass_name+"_"+width_name)

#print(histo_names)
os.system(f'root -l \'PlotVariablesMC.C(\"{allFiles}\", \"{level}\",\"{year}\", \"{histo_names}\", \"{mJJCut}\")\'')
    #break
