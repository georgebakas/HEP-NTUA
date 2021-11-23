import os, sys

years = ['2016_preVFP', '2016_postVFP', '2017', '2018']
suffix = ['out', 'err', 'log', 'sub']

for iy in years:
    for suf in suffix:
        os.system(f'rm {iy}/*.{suf}')