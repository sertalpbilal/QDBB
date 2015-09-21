
# This file gathers the results from test folder and put them into table

import os, os.path

result = []

# outer loop: for all files in system
directory = '../test/portfolio-output/'
files = os.listdir(directory)
for efile in sorted(files):
    ofile = open(directory + efile, 'r')
    found = 0
    exp = []
    for line in ofile:
        words = line.split()
        if "Optimal value:" in line:
            exp.append(words[3])
        if "Total time elapsed:" in line:
            exp.append(words[4])
        if "nodes generated:" in line:
            exp.append(words[5])
        if "nodes processed:" in line:
            exp.append(words[5])
        if "cuts applied:" in line:
            exp.append(words[4])
        if "SUMMARY" in line:
            found = 1
            exp.append(efile)
    if(found):
        result.append(exp)
    ofile.close()

for ex in result:
    print ex

#for a in range(1,20):
#    print '%02d' %(a)
