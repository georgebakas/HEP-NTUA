import glob
import os
import sys
import SubmitCondorJobs as scj
import subprocess


print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []
mJJCuts = [1000] #, 1200, 1400, 1600, 1800, 2000]
#output = subprocess.check_output(command_to_run,
#        shell=True, stdout=subprocess.PIPE)


command = f'xrdfs root://grid02.physics.uoi.gr ls -u /store/user/ipapakri/ttbar/MC/Signal/{year}'
output = subprocess.check_output(command, shell=True)

print(output)
print(type(output))
output = output.decode("utf-8")
split_files = output.split('\n')
print(split_files)

os.system('voms-proxy-init -voms cms --out /afs/cern.ch/user/g/gbakas/.globus/myProxy_certificate')

# deprecated for 2018, they have move to Ioannina
#for ifile, file_name in enumerate(glob.iglob('/eos/cms/store/user/ipapakri/ttbar/MC/Signal/{}/*.root'.format(year), recursive=True)):

for ifile, file_name in enumerate(split_files):
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

    file_name_to_send = file_name.split('/')[-1]

    for mjj_cut in mJJCuts:
        print('Current mjj cut:', mjj_cut)
        #os.system(f'root -l -b -q \'FillHistograms_Reduced_JES.C(\"{file_name}\",\"{split_file_underscore[0]}\", \"{jes_variation}\",\"{year}\",{mjj_cut})\'')
        #os.system(f'root -l -b -q \'FillHistograms_Extended_JES.C(\"{file_name}\",\"{split_file_underscore[0]}\", \"{jes_variation}\",\"{year}\",{mjj_cut})\'')
        argument = f'-l -b -q FillHistograms_Reduced_JES.C(\\\"{file_name_to_send}\\\",\\\"{split_file_underscore[0]}\\\",\\\"{jes_variation}\\\",\\\"{year}\\\",{mjj_cut}) {file_name}'
        output_file = f'HistoReduced_{mjj_cut}_{split_file_underscore[0]}_{jes_variation}.root'
        scj.submitCondorJobs('submit_.sh', argument, 'FillHistograms_Reduced_JES.C, TemplateConstants.h', output_file)

        argument = f'-l -b -q FillHistograms_Extended_JES.C(\\\"{file_name_to_send}\\\",\\\"{split_file_underscore[0]}\\\",\\\"{jes_variation}\\\",\\\"{year}\\\",{mjj_cut}) {file_name}'
        output_file_red = f'Histo_{mjj_cut}_{split_file_underscore[0]}_{jes_variation}.root'
        scj.submitCondorJobs('submit_.sh', argument, 'FillHistograms_Extended_JES.C, TemplateConstants.h', output_file_red)
        #break
    #break
