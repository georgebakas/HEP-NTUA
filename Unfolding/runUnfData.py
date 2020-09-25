import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
isParton = sys.argv[2]
iMethod = sys.argv[3]

os.system(f'root -l -b -q \'Unfold_data.cpp(\"{year}\", {isParton}, {iMethod})\'')
        #break
