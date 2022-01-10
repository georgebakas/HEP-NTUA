#!/bin/sh

# Response Matrices

# Nominal files

nohup python3 RunResponses_nominal.py 2016_preVFP > nom_responses_16pre.out &
nohup python3 RunResponses_nominal.py 2016_postVFP > nom_responses_16post.out &
nohup python3 RunResponses_nominal.py 2017 > nom_responses_17.out &
nohup python3 RunResponses_nominal.py 2018 > nom_responses_18.out &

# bTag Variation

nohup python3 RunResponses_bTag.py 2016_preVFP > btag_responses_16pre.out &
nohup python3 RunResponses_bTag.py 2016_postVFP > btag_responses_16post.out &
nohup python3 RunResponses_bTag.py 2017 > btag_responses_17.out &
nohup python3 RunResponses_bTag.py 2018 > btag_responses_18.out &

# PS Weights

nohup python3 RunResponses_PS_PDF.py 2016_preVFP PSWeights > ps_responses_16pre.out &
nohup python3 RunResponses_PS_PDF.py 2016_postVFP PSWeights > ps_responses_16post.out &
nohup python3 RunResponses_PS_PDF.py 2017 PSWeights > ps_responses_17.out &
nohup python3 RunResponses_PS_PDF.py 2018 PSWeights > ps_responses_18.out &

# PDF weights 

nohup python3 RunResponses_PS_PDF.py 2016_preVFP PDFWeights > pdf_responses_16pre.out &
nohup python3 RunResponses_PS_PDF.py 2016_postVFP PDFWeights > pdf_responses_16post.out &
nohup python3 RunResponses_PS_PDF.py 2017 PDFWeights > pdf_responses_17.out &
nohup python3 RunResponses_PS_PDF.py 2018 PDFWeights > pdf_responses_18.out &

# Scale Weights

nohup python3 RunResponses_PS_PDF.py 2016_preVFP ScaleWeights > scale_responses_16pre.out &
nohup python3 RunResponses_PS_PDF.py 2016_postVFP ScaleWeights > scale_responses_16post.out &
nohup python3 RunResponses_PS_PDF.py 2017 ScaleWeights > scale_responses_17.out &
nohup python3 RunResponses_PS_PDF.py 2018 ScaleWeights > scale_responses_18.out &

# JES Weights

#nohup python3 RunResponses_PS_PDF.py 2016_preVFP JES > jes_responses_16pre.out &
#nohup python3 RunResponses_PS_PDF.py 2016_postVFP JES > jes_responses_16post.out &
#nohup python3 RunResponses_PS_PDF.py 2017 JES > jes_responses_17.out &
#nohup python3 RunResponses_PS_PDF.py 2018 JES > jes_responses_18.out &