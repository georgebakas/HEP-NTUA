import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

for ifile, file_name in enumerate(glob.iglob('/eos/cms/store/user/gbakas/ZprimeToTT/mc-{}_btag/*.root'.format(year), recursive=True)):
    print('file: {}'.format(ifile))
    split_file_name = file_name.split('/')
    split_file_in_underscores = split_file_name[-1].split('_')
    mass_name = split_file_in_underscores[1]
    #print(mass_name)
    #print(split_file_name[-1])
    os.system(f'root -l -b -q \'FillHistograms_CutFlowJetMassSoftDrop.C(\"{split_file_name[-1]}\",\"{mass_name}\",\"{year}\")\'')
    #break
