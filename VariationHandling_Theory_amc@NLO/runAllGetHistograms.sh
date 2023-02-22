#!/bin/sh

# Response Matrices

# Nominal files

# nohup python3 RunGetHistograms_nominal.py 2016_preVFP > nom_histograms_16pre.out &
# nohup python3 RunGetHistograms_nominal.py 2016_postVFP > nom_histograms_16post.out &
# nohup python3 RunGetHistograms_nominal.py 2017 > nom_histograms_17.out &
nohup python3 RunGetHistograms_nominal.py 2018 > nom_histograms_18.out &

# bTag Variation

#nohup python3 RunGetHistograms_bTag.py 2016_preVFP > btag_histograms_16pre.out &
#nohup python3 RunGetHistograms_bTag.py 2016_postVFP > btag_histograms_16post.out &
#nohup python3 RunGetHistograms_bTag.py 2017 > btag_histograms_17.out &
#nohup python3 RunGetHistograms_bTag.py 2018 > btag_histograms_18.out &

# PS Weights
# No event counter ps weights
#nohup python3 RunGetHistograms_PS_PDF.py 2016_preVFP PSWeights > ps_histograms_16pre.out &
#nohup python3 RunGetHistograms_PS_PDF.py 2016_postVFP PSWeights > ps_histograms_16post.out &
#nohup python3 RunGetHistograms_PS_PDF.py 2017 PSWeights > ps_histograms_17.out &
nohup python3 RunGetHistograms_PS_PDF.py 2018 PSWeights > ps_histograms_18.out &

# PDF weights 

#nohup python3 RunGetHistograms_PS_PDF.py 2016_preVFP PDFWeights > pdf_histograms_16pre.out &
#nohup python3 RunGetHistograms_PS_PDF.py 2016_postVFP PDFWeights > pdf_histograms_16post.out &
#nohup python3 RunGetHistograms_PS_PDF.py 2017 PDFWeights > pdf_histograms_17.out &
nohup python3 RunGetHistograms_PS_PDF.py 2018 PDFWeights > pdf_histograms_18.out &

# Scale Weights

#nohup python3 RunGetHistograms_PS_PDF.py 2016_preVFP ScaleWeights > scale_histograms_16pre.out &
#nohup python3 RunGetHistograms_PS_PDF.py 2016_postVFP ScaleWeights > scale_histograms_16post.out &
#nohup python3 RunGetHistograms_PS_PDF.py 2017 ScaleWeights > scale_histograms_17.out &
nohup python3 RunGetHistograms_PS_PDF.py 2018 ScaleWeights > scale_histograms_18.out &
