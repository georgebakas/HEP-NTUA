import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]

allFiles = []
histo_names = []

for ifile, file_name in enumerate(glob.iglob('{}/HistoCutFlowJetMassSoftDrop_*.root'.format(year), recursive=True)):
    print('file: {}'.format(file_name))
    file_to_append = file_name.split('/')
    split_file_in_underscores = file_to_append[-1].split('_')
    mass_name = split_file_in_underscores[2]
    width_name = split_file_in_underscores[3]
    histo_names = mass_name+"_"+width_name
    #os.system(f'root -l -q \'PlotCutFlow.C(\"{year}\", \"{file_to_append[-1]}\",\"{histo_names}\")\'')
    os.system(f'root -l -b -q \'PlotCutFlow_JetMassSoftDrop.C(\"{year}\", \"{file_to_append[-1]}\",\"{histo_names}\")\'')
    #break
