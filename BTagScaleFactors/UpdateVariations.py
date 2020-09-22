import glob 
import os
import sys 

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

for ifile, file_name in enumerate(glob.iglob('/eos/cms/store/user/ipapakri/ttbar/MC/Signal/{}/variations/*.root'.format(year), recursive=True)):
	#print(file_name)
	#print(f'root -l -b -q \'myMacro.C(\"{file_name}\")\'')
	if ifile > 21:
		print("file {}".format(ifile))
		split_file_name = file_name.split('/')
		print(split_file_name[-1])
		os.system(f'root -l -b -q \'UpdateBTagSF_Variations.C(\"{split_file_name[-1]}\", \"{year}\")\'')
	#break
