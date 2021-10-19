import glob
import os
import sys
import SubmitCondorJobs as scj
import subprocess

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []
mJJCuts = [1000]#, 1200, 1400, 1600, 1800, 2000]

command = f'xrdfs root://grid02.physics.uoi.gr ls -u /store/user/ipapakri/ttbar/MC/Signal/{year}/'
output = subprocess.check_output(command, shell=True)

#print(output)
#print(type(output))
output = output.decode("utf-8")
split_files = output.split('\n')
print(len(split_files))
os.system('voms-proxy-init -voms cms --out /afs/cern.ch/user/g/gbakas/.globus/myProxy_certificate')

for ifile, file_name in enumerate(split_files):
    if ('mtop' or 'hdamp') in file_name:
        print('file: {}'.format(file_name))
        split_file_name = file_name.split('/')
        split_file_underscore = split_file_name[-1].split('_')
        
        print(split_file_name[-1]) #this is the file name without the '/'

        print(split_file_underscore[0])
        # get the JES variation name
        mass_variation = split_file_underscore[1]
        
            
        print(mass_variation)
        
        file_name_to_send = file_name.split('/')[-1]

        for mjj_cut in mJJCuts:

            print('Current mjj cut:', mjj_cut)
            argument = f'-l -b -q ResponseMatrices_SystematicsFiles.C(\\\"{file_name_to_send}\\\",\\\"{split_file_underscore[0]}\\\",\\\"{mass_variation}\\\",\\\"{year}\\\",{mjj_cut}) {file_name}'
            print(argument)
            #               ResponsesEfficiency_TTTo2L2Nu_mtop166p5
            output_file = f'ResponsesEfficiency_{split_file_underscore[0]}_{mass_variation}.root'
            print(output_file)
            scj.submitCondorJobs('submit_.sh', argument, 'ResponseMatrices_SystematicsFiles.C, TemplateConstants.h', output_file)



