import sys
import os 

phase_space = ['Parton', 'Particle']

for p in phase_space:

    # Chi2 results file 
    norm = 'false'
    os.system(f'root -l -b -q \'ChiSquare.cpp(\"\", \"{p}\")\'')
    
    # normalised 
    os.system(f'root -l -b -q \'ChiSquare_normalised.cpp(\"\", \"{p}\")\'')