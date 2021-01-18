import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
allFiles = []

masses = [2000, 2500, 3000, 3500, 4000]
width = [0.01, 0.1, 0.3]

for imass in masses:
    if year == '2016' and imass == 3500:
        continue

    for w in width:
        iwidth = int(imass * w)

        if year == '2016' and imass == 2500 and w == 0.3:
            continue
        print('mass', imass, 'width:', iwidth)
        os.system(f'root -l -b -q \'plotStackSensitivity.C(\"{year}\",1600,{imass},{iwidth})\'')
