import glob
import os
import sys
import SubmitCondorJobs as scj

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []
mJJCuts = [1000]

for ifile, file_name in enumerate(glob.iglob('/eos/cms/store/user/ipapakri/ttbar/MC/Signal/{}/*.root'.format(year), recursive=True)):
    print('file: {}'.format(file_name))
    split_file_name = file_name.split('/')
    split_file_underscore = split_file_name[-1].split('_')
    print(split_file_underscore[0])
    print(split_file_name[-1])
    print(split_file_underscore[-1])
    # get the JES variation name
    if len(split_file_underscore) > 3:
        jes_variation = split_file_underscore[-1].split('.')[0]
    else:
        continue
    print(jes_variation)

    for mjj_cut in mJJCuts:
        print('Current mjj cut:', mjj_cut)
        #os.system(f'root -l -b -q \'FillHistograms_Reduced_JES.C(\"{file_name}\",\"{split_file_underscore[0]}\", \"{jes_variation}\",\"{year}\",{mjj_cut})\'')
        #os.system(f'root -l -b -q \'FillHistograms_Extended_JES.C(\"{file_name}\",\"{split_file_underscore[0]}\", \"{jes_variation}\",\"{year}\",{mjj_cut})\'')
        argument = f'-l -b -q ResponsesMatrices_JES.C(\\\"{file_name}\\\",\\\"{split_file_underscore[0]}\\\",\\\"{jes_variation}\\\",\\\"{year}\\\",{mjj_cut})'
        output_file = f'ResponsesEfficiency_{split_file_underscore[0]}_{jes_variation}.root'
        scj.submitCondorJobs('submit_.sh', argument, 'ResponsesMatrices_JES.C, TemplateConstants.h', output_file)
        #break
    #break
