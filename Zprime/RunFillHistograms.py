import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

xsec16 = {
 'M-1000':0.5392, 'M-1500':0.07238, 'M-2000':0.01375, 'M-2500':0.003, 'M-3000':0.0007569, 'M-3500':0.0001979, 'M-4000':5.484*pow(10,-5), 'M-5000':5.043*pow(10,-6)
}
xsec17 = {
 'M-1000':0.5392, 'M-1500':0.07238, 'M-2000':0.01375, 'M-2500':0.003, 'M-3000':0.0007569, 'M-3500':0.0001979, 'M-4000':5.484*pow(10,-5), 'M-5000':5.043*pow(10,-6)
}
xsec18 = {
 'M-1000':0.5392, 'M-1500':0.07238, 'M-2000':0.01375, 'M-2500':0.003, 'M-3000':0.0007569, 'M-3500':0.0001979, 'M-4000':5.484*pow(10,-5), 'M-5000':5.043*pow(10,-6)
}



for ifile, file_name in enumerate(glob.iglob('/eos/cms/store/user/gbakas/ZprimeToTT/mc-{}_btag/*.root'.format(year), recursive=True)):
	#print(file_name)
	print("file {}".format(ifile))
	split_file_name = file_name.split('/')
    split_file_in_underscores = file_name.split('_')
    mass_name = split_file_in_underscores[1]
    print(mass_name)

	print(split_file_name[-1])
	#print(file_name)
	#os.system(f'root -l -b -q \'UpdateBTagSF.C(\"{split_file_name[-1]}\", \"{year}\")\'')
	break
