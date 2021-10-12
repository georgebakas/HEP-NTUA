import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')


variations = ['Nominal', 'bTagVariation', 'JES', 'PDFWeights', 'PSWeights', 'ScaleWeights']
type_ = ['Parton', 'Particle']

for var in variations:
    for type_i in type_:
        os.system(f'root -l -b -q \'CombineEfficiency.cxx(\"{var}\",\"{type_i}\")\'')
        os.system(f'root -l -b -q \'CombineAcceptance.cxx(\"{var}\",\"{type_i}\")\'')
        #break
    #break
