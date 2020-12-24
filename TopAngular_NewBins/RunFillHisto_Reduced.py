import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

mass_ranges = [1000, 1200, 1400, 1600, 1800, 2000]

for imass in mass_ranges:
    if imass == 1000:
        continue
    print('mTTbar cut running now is: ', imass, ' GeV')
    for isel in range(0,5):
        if isel == 1:
            continue
        print('Selection is: ', isel)
        print('year is: ', year)
        os.system(f'root -l -b -q \'FillHistograms_Reduced.C(\"{year}\",{isel},{imass})\'')
