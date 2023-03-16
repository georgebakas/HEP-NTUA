import glob
import os
import sys

allFiles = []

ps_weights = {2:"isrRedHi", 3:"fsrRedHi", 4:"isrRedLo", 5:"fsrRedLo", 6:"isrDefHi", 7:"fsrDefHi",
            8:"isrDefLo", 9:"fsrDefLo", 10:"isrConHi", 11:"fsrConHi", 12:"isrConLo", 13:"fsrConLo"}

def unfold(year, weightType, isParton):
    
    for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/*TTToHadronic*.root'.format(year,weightType), recursive=True)):
        print(file_name)
        print(weightType)
        if weightType != 'SystematicsFiles':
            split_file_name = file_name.split('/')
            split_file_underscore = split_file_name[-1].split('_')

            print(split_file_underscore[1])
            tt_process_name = split_file_underscore[1]

            dot_split = split_file_underscore[-1].split('.')
            weight_sufix = dot_split[0]
            print('weight_suffix', weight_sufix)

            #os.system(f'root -l -b -q \'Unfold_data.cpp(\"{year}\",\"{weightType}\",\"{weight_sufix}\", true)\'')
        elif weightType == 'Nominal':
            weight_sufix = ""
        else:
            #TT, and TTJets files are handled below
            split_file_name = file_name.split('/')
            split_file_underscore = split_file_name[-1].split('_')

            print(split_file_underscore[1])
            tt_process_name = split_file_underscore[1]

            dot_split = split_file_name[-1].replace(split_file_underscore[0]+'_TTToHadronic_', '')
            dot_split = dot_split.split('.')
            weight_sufix = dot_split[0]
            print('weight_suffix', weight_sufix)
            #(TString inYear, TString dir, TString inputFile, bool isThreeProcesses, bool isParton = true, int unfoldMethod=1)

        os.system(f'root -l -b -q \'Unfold_Combined.cxx(\"{weightType}\",\"{weight_sufix}\",\"{isParton}\")\'')
        #break

    '''
    combined_files = ['TT','TTJets']
    #check for files of type TT_Tune
    if weightType == 'SystematicsFiles':
        for tt in combined_files:
            for ifile, file_name in enumerate(glob.iglob('{}/Responses{}/*{}_*.root'.format(year,weightType, tt), recursive=True)):
                split_file_name = file_name.split('/')
                split_file_underscore = split_file_name[-1].split('_')

                print(split_file_underscore[1])
                tt_process_name = split_file_underscore[1]

                dot_split = split_file_name[-1].replace('ResponsesEfficiency_', '')
                dot_split = dot_split.split('.')
                weight_sufix = dot_split[0]
                print('weight_suffix', weight_sufix)
                os.system(f'root -l -b -q \'Unfold_Combined.C(\"{year}\",\"{weightType}\",\"{weight_sufix}\", false)\'')
    #eof
    '''

if __name__ == '__main__':
    
    types = ['Nominal', 'NominalTopSF','bTagVariation', 'topTaggingVariation', 'PSWeights', 'ScaleWeights', 'PDFWeights', 'JES']
    types = ['NominalTopSF']
    year = '2017'
    parton_particle_types = ['Parton', 'Particle']
    
    #unfold('2017', 'PDFWeights', 'true')
    #exit()
    for itype in types:
        print(itype)
        #if 'ScaleWeights' in itype: 
        unfold(year, itype, 'true')
        unfold(year, itype, 'false')
        #break