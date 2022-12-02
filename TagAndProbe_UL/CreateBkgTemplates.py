import glob
import os
import sys

print(f'Number of arguments {len(sys.argv)}')
print(f'Argument list: {str(sys.argv)}')

year = sys.argv[1]
tag_probe = ["probe", "SR"]

for tg in tag_probe:
    os.system(f'root -l -q -b \'CreateBkgTemplates.C(\"{year}\", \"{tg}\")\'')
    #break
