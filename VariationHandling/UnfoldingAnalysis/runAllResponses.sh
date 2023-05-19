#!/bin/sh

# # Response Matrices

# # Nominal files

nohup python3 RunResponses_nominal.py 2016_preVFP > nom_responses_16pre.out &
nohup python3 RunResponses_nominal.py 2016_postVFP > nom_responses_16post.out &
nohup python3 RunResponses_nominal.py 2017 > nom_responses_17.out &
nohup python3 RunResponses_nominal.py 2018 > nom_responses_18.out &
