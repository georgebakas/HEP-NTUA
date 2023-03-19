#!/bin/sh

# # Extended files needed for the analysis

# # Nominal files

# nohup python3 RunFillHisto_nominal.py 2016_preVFP True > nom_ext_16pre.out &
# nohup python3 RunFillHisto_nominal.py 2016_postVFP True > nom_ext_16post.out &
# nohup python3 RunFillHisto_nominal.py 2017 True > nom_ext_17.out &
# nohup python3 RunFillHisto_nominal.py 2018 True > nom_ext_18.out &

# Nominal Files with Top Tagger SF (Extended)

nohup python3 RunFillHisto_nominalTopSF.py 2016_preVFP True > NomTopTag_ext_16pre.out &
nohup python3 RunFillHisto_nominalTopSF.py 2016_postVFP True > NomTopTag_ext_16post.out &
nohup python3 RunFillHisto_nominalTopSF.py 2017 True > NomTopTag_ext_17.out &
nohup python3 RunFillHisto_nominalTopSF.py 2018 True > NomTopTag_ext_18.out &

# # bTag Variation

# nohup python3 RunFillHisto_bTag.py 2016_preVFP True> btag_ext_16pre.out &
# nohup python3 RunFillHisto_bTag.py 2016_postVFP True > btag_ext_16post.out &
# nohup python3 RunFillHisto_bTag.py 2017 True > btag_ext_17.out &
# nohup python3 RunFillHisto_bTag.py 2018 True > btag_ext_18.out &

# # TopTagger SF Variation

nohup python3 RunFillHisto_topTag.py 2016_preVFP True> topTag_ext_16pre.out &
nohup python3 RunFillHisto_topTag.py 2016_postVFP True > topTag_ext_16post.out &
nohup python3 RunFillHisto_topTag.py 2017 True > topTag_ext_17.out &
nohup python3 RunFillHisto_topTag.py 2018 True > topTag_ext_18.out &

# # PS Weights

# nohup python3 RunFillHisto_PS_PDF.py 2016_preVFP PSWeights True> ps_ext_16pre.out &
# nohup python3 RunFillHisto_PS_PDF.py 2016_postVFP PSWeights True > ps_ext_16post.out &
# nohup python3 RunFillHisto_PS_PDF.py 2017 PSWeights True > ps_ext_17.out &
# nohup python3 RunFillHisto_PS_PDF.py 2018 PSWeights True > ps_ext_18.out &

# # PDF weights 

# nohup python3 RunFillHisto_PS_PDF.py 2016_preVFP PDFWeights True> pdf_ext_16pre.out &
# nohup python3 RunFillHisto_PS_PDF.py 2016_postVFP PDFWeights True > pdf_ext_16post.out &
# nohup python3 RunFillHisto_PS_PDF.py 2017 PDFWeights True > pdf_ext_17.out &
# nohup python3 RunFillHisto_PS_PDF.py 2018 PDFWeights True > pdf_ext_18.out &

# # Scale Weights

# nohup python3 RunFillHisto_PS_PDF.py 2016_preVFP ScaleWeights True > scale_ext_16pre.out &
# nohup python3 RunFillHisto_PS_PDF.py 2016_postVFP ScaleWeights True > scale_ext_16post.out &
# nohup python3 RunFillHisto_PS_PDF.py 2017 ScaleWeights True > scale_ext_17.out &
# nohup python3 RunFillHisto_PS_PDF.py 2018 ScaleWeights True > scale_ext_18.out &

# # JES Weights

# #nohup python3 RunFillHisto_JES.py 2016_preVFP JES True > jes_ext_16pre.out &
# #nohup python3 RunFillHisto_JES.py 2016_postVFP JES True > jes_ext_16post.out &
# #nohup python3 RunFillHisto_JES.py 2017 JES True > jes_ext_17.out &
# #nohup python3 RunFillHisto_JES.py 2018 JES True > jes_ext_18.out &

# # Mass Variation Files 

# #nohup python3 RunFillHisto_JES.py 2016_preVFP JES True > scale_ext_16pre.out &
# #nohup python3 RunFillHisto_JES.py 2016_postVFP JES True > scale_ext_16post.out &
# #nohup python3 RunFillHisto_JES.py 2017 JES True > scale_ext_17.out &
# #nohup python3 RunFillHisto_JES.py 2018 JES True > scale_ext_18.out &

# # Reduced files needed for the analysis

# # Nominal files

# nohup python3 RunFillHisto_nominal.py 2016_preVFP False > nom_red_16pre.out &
# nohup python3 RunFillHisto_nominal.py 2016_postVFP False > nom_red_16post.out &
# nohup python3 RunFillHisto_nominal.py 2017 False > nom_red_17.out &
# nohup python3 RunFillHisto_nominal.py 2018 False > nom_red_18.out &

# Nominal Files with Top Tagger SF (Extended)

nohup python3 RunFillHisto_nominalTopSF.py 2016_preVFP False> NomTopTag_red_16pre.out &
nohup python3 RunFillHisto_nominalTopSF.py 2016_postVFP False > NomTopTag_red_16post.out &
nohup python3 RunFillHisto_nominalTopSF.py 2017 False > NomTopTag_red_17.out &
nohup python3 RunFillHisto_nominalTopSF.py 2018 False > NomTopTag_red_18.out &

# # bTag Variation

# nohup python3 RunFillHisto_bTag.py 2016_preVFP False> btag_red_16pre.out &
# nohup python3 RunFillHisto_bTag.py 2016_postVFP False > btag_red_16post.out &
# nohup python3 RunFillHisto_bTag.py 2017 False > btag_red_17.out &
# nohup python3 RunFillHisto_bTag.py 2018 False > btag_red_18.out &

# # TopTagger SF Variation

nohup python3 RunFillHisto_topTag.py 2016_preVFP False > topTag_red_16pre.out &
nohup python3 RunFillHisto_topTag.py 2016_postVFP False > topTag_red_16post.out &
nohup python3 RunFillHisto_topTag.py 2017 False > topTag_red_17.out &
nohup python3 RunFillHisto_topTag.py 2018 False > topTag_red_18.out &

# # PS Weights

# nohup python3 RunFillHisto_PS_PDF.py 2016_preVFP PSWeights False> ps_red_16pre.out &
# nohup python3 RunFillHisto_PS_PDF.py 2016_postVFP PSWeights False > ps_red_16post.out &
# nohup python3 RunFillHisto_PS_PDF.py 2017 PSWeights False > ps_red_17.out &
# nohup python3 RunFillHisto_PS_PDF.py 2018 PSWeights False > ps_red_18.out &

# # PDF weights 

# nohup python3 RunFillHisto_PS_PDF.py 2016_preVFP PDFWeights False> pdf_red_16pre.out &
# nohup python3 RunFillHisto_PS_PDF.py 2016_postVFP PDFWeights False > pdf_red_16post.out &
# nohup python3 RunFillHisto_PS_PDF.py 2017 PDFWeights False > pdf_red_17.out &
# nohup python3 RunFillHisto_PS_PDF.py 2018 PDFWeights False > pdf_red_18.out &

# # Scale Weights

# nohup python3 RunFillHisto_PS_PDF.py 2016_preVFP ScaleWeights False > scale_red_16pre.out &
# nohup python3 RunFillHisto_PS_PDF.py 2016_postVFP ScaleWeights False > scale_red_16post.out &
# nohup python3 RunFillHisto_PS_PDF.py 2017 ScaleWeights False > scale_red_17.out &
# nohup python3 RunFillHisto_PS_PDF.py 2018 ScaleWeights False > scale_red_18.out &

# # JES Weights

# #nohup python3 RunFillHisto_JES.py 2016_preVFP JES False > jes_red_16pre.out &
# #nohup python3 RunFillHisto_JES.py 2016_postVFP JES False > jes_red_16post.out &
# #nohup python3 RunFillHisto_JES.py 2017 JES False > jes_red_17.out &
# #nohup python3 RunFillHisto_JES.py 2018 JES False > jes_red_18.out &

# # Mass Variation Files

# #nohup python3 RunFillHisto_JES.py 2016_preVFP JES False > scale_red_16pre.out &
# #nohup python3 RunFillHisto_JES.py 2016_postVFP JES False > scale_red_16post.out &
# #nohup python3 RunFillHisto_JES.py 2017 JES False > scale_red_17.out &
# #nohup python3 RunFillHisto_JES.py 2018 JES False > scale_red_18.out &