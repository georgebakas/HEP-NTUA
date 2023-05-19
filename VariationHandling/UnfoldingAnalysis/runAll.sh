#!/bin/sh

# # Extended files needed for the analysis

# # Nominal files

nohup python3 RunFillHisto_nominal.py 2016_preVFP True > nom_ext_16pre.out &
nohup python3 RunFillHisto_nominal.py 2016_postVFP True > nom_ext_16post.out &
nohup python3 RunFillHisto_nominal.py 2017 True > nom_ext_17.out &
nohup python3 RunFillHisto_nominal.py 2018 True > nom_ext_18.out &

# # Reduced files needed for the analysis

# # Nominal files

nohup python3 RunFillHisto_nominal.py 2016_preVFP False > nom_red_16pre.out &
nohup python3 RunFillHisto_nominal.py 2016_postVFP False > nom_red_16post.out &
nohup python3 RunFillHisto_nominal.py 2017 False > nom_red_17.out &
nohup python3 RunFillHisto_nominal.py 2018 False > nom_red_18.out &