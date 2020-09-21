import glob 
import os
import sys 

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]

for file_name in glob.iglob('/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2016/variations/*.root', recursive=True):
	#print(file_name)
	#print(f'root -l -b -q \'myMacro.C(\"{file_name}\")\'')
	split_file_name = file_name.split('/')
	print(split_file_name[-1])
	os.system(f'root -l -b -q \'UpdateBTagSF_Variations.C(\"{split_file_name[-1]}\", \"{year}\")\'')
	break
