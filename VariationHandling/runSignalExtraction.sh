#!/bin/sh

# Create Mass Fit Results

# Nominal files

python3 SignalExtraction.py 2016_preVFP Nominal
python3 SignalExtraction.py 2016_postVFP Nominal
python3 SignalExtraction.py 2017 Nominal
python3 SignalExtraction.py 2018 Nominal

# bTagVariation

python3 SignalExtraction.py 2016_preVFP bTagVariation
python3 SignalExtraction.py 2016_postVFP bTagVariation
python3 SignalExtraction.py 2017 bTagVariation
python3 SignalExtraction.py 2018 bTagVariation

# topTaggingVariation

python3 SignalExtraction.py 2016_preVFP topTaggingVariation
python3 SignalExtraction.py 2016_postVFP topTaggingVariation
python3 SignalExtraction.py 2017 topTaggingVariation
python3 SignalExtraction.py 2018 topTaggingVariation

# PS Weights

python3 SignalExtraction.py 2016_preVFP PSWeights
python3 SignalExtraction.py 2016_postVFP PSWeights
python3 SignalExtraction.py 2017 PSWeights
python3 SignalExtraction.py 2018 PSWeights

# Scale Weights 

python3 SignalExtraction.py 2016_preVFP ScaleWeights
python3 SignalExtraction.py 2016_postVFP ScaleWeights
python3 SignalExtraction.py 2017 ScaleWeights
python3 SignalExtraction.py 2018 ScaleWeights

# PDF Weights

python3 SignalExtraction.py 2016_preVFP PDFWeights
python3 SignalExtraction.py 2016_postVFP PDFWeights
python3 SignalExtraction.py 2017 PDFWeights
python3 SignalExtraction.py 2018 PDFWeights


# JES Weights

python3 SignalExtraction.py 2016_preVFP JES
python3 SignalExtraction.py 2016_postVFP JES
python3 SignalExtraction.py 2017 JES
python3 SignalExtraction.py 2018 JES

# Mass Variation and hdamp files 

python3 SignalExtraction.py 2016_preVFP SystematicsFiles
python3 SignalExtraction.py 2016_postVFP SystematicsFiles
python3 SignalExtraction.py 2017 SystematicsFiles
python3 SignalExtraction.py 2018 SystematicsFiles