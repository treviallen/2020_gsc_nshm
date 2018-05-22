
from numpy import array, savetxt, vstack, hstack
from sys import argv
from os import path

table1 = argv[1]
table2 = argv[2]

# parse tables
lines1 = open(table1).readlines()
lines2 = open(table2).readlines()

# make header
header = ' '.join(('Difference between', path.split(table2)[-1], 'and', path.split(table1)[-1], 'ground motions in regular (non-log) units')) + '\n'
header += lines1[1]
header += lines1[2]
header += lines1[3].strip('\n')

i = 0
for line1, line2 in zip(lines1[4:], lines2[4:]):
    dat1str = line1.strip().split()
    magdist = array([float(val) for val in dat1str[0:2]])
    	
    dat1 = array([float(val) for val in dat1str[2:]])
    
    dat2str = line2.strip().split()
    dat2 = array([float(val) for val in dat2str[2:]])
    
    diff = 10**(dat2 - dat1)
    
    newline = hstack((magdist, diff))
    	
    if i == 0:
        diffArray = newline
    else:
        diffArray = vstack((diffArray, newline))
        
    i += 1
       
# make file name
part1 = path.split(table1)[-1].split('.')
part2 = path.split(table2)[-1].split('.')

newFile = '.'.join((part1[0], part1[1], part2[1], part1[2],'txt'))
outFolder = 'gmm_amp_diff'
newPath = path.join(outFolder, newFile)

# write file
savetxt(newPath, diffArray, fmt='%.3f', delimiter='\t', header=header)

