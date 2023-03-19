import os

types = ['Nominal', 'NominalTopSF', 'bTagVariation', 'topTaggingVariation', 'PSWeights', 'ScaleWeights', 'PDFWeights', 'JES']
types = ['NominalTopSF', 'topTaggingVariation']

parton_particle_types = ['Parton', 'Particle']

for weightType in types:
    for parton_particle in parton_particle_types:
        print('--------------------------')
        print(weightType, parton_particle)
        os.system(f'root -l -b -q \'CombineEfficiency.cxx(\"{weightType}\",\"{parton_particle}\")\'')
        os.system(f'root -l -b -q \'CombineAcceptance.cxx(\"{weightType}\",\"{parton_particle}\")\'')
        #break
