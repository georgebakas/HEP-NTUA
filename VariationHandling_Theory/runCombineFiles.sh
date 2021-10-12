#!/bin/sh

# Create Combined Files

# Nominal files

python3 combineFiles.py 2016_preVFP Nominal True
python3 combineFiles.py 2016_postVFP Nominal True
python3 combineFiles.py 2017 Nominal True
python3 combineFiles.py 2018 Nominal True

# bTagVariation

python3 combineFiles.py 2016_preVFP bTagVariation True
python3 combineFiles.py 2016_postVFP bTagVariation True
python3 combineFiles.py 2017 bTagVariation True
python3 combineFiles.py 2018 bTagVariation True

# PS Weights

python3 combineFiles.py 2016_preVFP PSWeights True
python3 combineFiles.py 2016_postVFP PSWeights True
python3 combineFiles.py 2017 PSWeights True
python3 combineFiles.py 2018 PSWeights True

# Scale Weights 

python3 combineFiles.py 2016_preVFP ScaleWeights True
python3 combineFiles.py 2016_postVFP ScaleWeights True
python3 combineFiles.py 2017 ScaleWeights True
python3 combineFiles.py 2018 ScaleWeights True

# PDF Weights

python3 combineFiles.py 2016_preVFP PDFWeights True
python3 combineFiles.py 2016_postVFP PDFWeights True
python3 combineFiles.py 2017 PDFWeights True
python3 combineFiles.py 2018 PDFWeights True

# Mass Variation Files 

python3 combineFiles.py 2016_preVFP SystematicsFiles True
python3 combineFiles.py 2016_postVFP SystematicsFiles True
python3 combineFiles.py 2017 SystematicsFiles True
python3 combineFiles.py 2018 SystematicsFiles True


#!/bin/sh

# Create Combined Files (Reduced)

# Nominal files

python3 combineFiles.py 2016_preVFP Nominal False
python3 combineFiles.py 2016_postVFP Nominal False
python3 combineFiles.py 2017 Nominal False
python3 combineFiles.py 2018 Nominal False

# bTagVariation

python3 combineFiles.py 2016_preVFP bTagVariation False
python3 combineFiles.py 2016_postVFP bTagVariation False
python3 combineFiles.py 2017 bTagVariation False
python3 combineFiles.py 2018 bTagVariation False

# PS Weights

python3 combineFiles.py 2016_preVFP PSWeights False
python3 combineFiles.py 2016_postVFP PSWeights False
python3 combineFiles.py 2017 PSWeights False
python3 combineFiles.py 2018 PSWeights False

# EXTENDED 

# Scale Weights 

python3 combineFiles.py 2016_preVFP ScaleWeights False
python3 combineFiles.py 2016_postVFP ScaleWeights False
python3 combineFiles.py 2017 ScaleWeights False
python3 combineFiles.py 2018 ScaleWeights False

# PDF Weights

python3 combineFiles.py 2016_preVFP PDFWeights False
python3 combineFiles.py 2016_postVFP PDFWeights False
python3 combineFiles.py 2017 PDFWeights False
python3 combineFiles.py 2018 PDFWeights False

# JES Weights

python3 combineFiles.py 2016_preVFP JES False
python3 combineFiles.py 2016_postVFP JES False
python3 combineFiles.py 2017 JES False
python3 combineFiles.py 2018 JES False

# Mass Variation Files 

python3 combineFiles.py 2016_preVFP SystematicsFiles False
python3 combineFiles.py 2016_postVFP SystematicsFiles False
python3 combineFiles.py 2017 SystematicsFiles False
python3 combineFiles.py 2018 SystematicsFiles False