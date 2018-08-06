import sys

file1 = sys.argv[1]
outfile = open(file1.split('.')[0]+'.tab.txt','w')

with open(file1) as F:
	for line in F:
		line = line.strip('\n').split()
		outfile.write('\t'.join(line)+'\n')
