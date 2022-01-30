import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

masses = [1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500]
width = [0.01]#, 0.1, 0.3]

for imass in masses:

    for w in width:
        iwidth = int(imass * w)

        if year == '2016' and imass == 2500 and w == 0.3:
            continue
        print('mass', imass, 'width:', iwidth)

        os.system(f'root -l -b -q \'plotStackHisto.C(\"{year}\",1000,{imass},{iwidth})\'')
        os.system(f'root -l -b -q \'plotStackHisto_ZprimeInc.C(\"{year}\",1000,{imass},{iwidth})\'')
        #os.system(f'root -l -b -q \'plotStackHisto_DeltaY.C(\"{year}\",1400,{imass},{iwidth})\'')
