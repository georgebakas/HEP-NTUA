import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

for ifile, file_name in enumerate(glob.iglob('/eos/cms/store/user/ipapakri/ttbar/MC/Signal/{}/variations/TTTo*.root'.format(year), recursive=True)):
	#print(file_name)
	print("file {}".format(ifile))
	split_file_name = file_name.split('/')
	print(split_file_name[-1])
	#print(file_name)
	os.system(f'root -l -b -q \'UpdateBTagSF.C(\"{split_file_name[-1]}\", \"{year}\")\'')
	#break
