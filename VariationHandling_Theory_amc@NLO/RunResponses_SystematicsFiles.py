import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

for ifile, file_name in enumerate(glob.iglob('/eos/cms/store/user/gbakas/ttbar/topTagger/ul-{}/variations/combined/*.root'.format(year), recursive=True)):
    print('file: {}'.format(ifile))
    split_file_name = file_name.split('/')
    split_file_underscore = split_file_name[-1].split('_')
    print('process:', split_file_underscore[0])
    print(split_file_name[-1])
    os.system(f'root -l -b -q \'ResponseMatrices_SystematicsFiles.C(\"{split_file_name[-1]}\",\"{split_file_underscore[0]}\",\"{year}\")\'')
    #break
