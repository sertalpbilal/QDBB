
# This create.py script search the data folder and
# create condor submission file (condor.sub) for all data

import os

# Step -1: Remove all files in test portfolio-output folder

directory = '../test/portfolio-output/'
files = os.listdir(directory)
for efile in files:
    os.remove(directory + efile)



# Step 0: Condor for AA problems

# Step 0.0: Arrays
capital = [10000]
ret = [0.01, 0.04]
cardinality = [5, 10]
branch = ['mf','hc'] # random is out
cut = ['mf', 'hc'] # random is out
search = ['df0', 'df1', 'breadth', 'best']
cutiter = [1]
cutperiter = [1]

# Step 0.1: Open file
cfile = open('condor.sub','w')
ifile = open('exp.txt','w')

# Step 0.2: Print common values
common_command = \
'Executable = ../test/portfolio \n\
Universe   = vanilla\n\
getenv     = true\n\
transfer_executable = false \n\n'
cfile.write(common_command)

# Step 0.3: Loop over arrays

comb = [(c,r,br,cu,sea) for c in capital for r in ret for br in branch for cu in cut for sea in search]
index = 0

for (c,r,br,cu,sea) in comb:
    index = index + 1
    ci = 0
    cp = 0
    # run command
    run_command =  \
'arguments  = -type roundlot -d AA -a 20 -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -C %d -r %f -f 0 -o 2\n\
output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, c, r, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d - round\tb:%s c:%s s:%s i:%d p:%d C:%d r:%.2f\n' %(index, br, cu, sea, ci, cp, c, r)
    ifile.write(prob_index)

    # x = 1
    index += 1
    run_command =  \
'arguments  = -type roundlot -d AA -a 20 -b %s -c %s -s %s -x 1 -l 100 -i %d -p %d -C %d -r %f -f 0 -o 2\n\
output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, c, r, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d - round\tb:%s c:%s s:%s i:%d p:%d C:%d r:%.2f x:1\n' %(index, br, cu, sea, ci, cp, c, r)
    ifile.write(prob_index)
    
    comb2 = [(ci,cp) for ci in cutiter for cp in cutperiter]
    for (ci,cp) in comb2:
        # run command
        index += 1
        run_command =  \
'arguments  = -type roundlot -d AA -a 20 -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -C %d -r %f -f 0 -o 2\n\
output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, c, r, index)
        # write to file
        cfile.write(run_command)
        prob_index = '%04d - round\tb:%s c:%s s:%s i:%d p:%d C:%d r:%.2f\n' %(index, br, cu, sea, ci, cp, c, r)
        ifile.write(prob_index)

#  -----------  CARDINALITY -----------

comb3 = [(r,br,cu,sea,k) for r in ret for br in branch for cu in cut for sea in search for k in cardinality]
for (r,br,cu,sea,k) in comb3:
    index = index+1
    ci = 0
    cp = 0
    run_command =  \
    'arguments  = -type cardinality -d AA -a 20 -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
    output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, r, k, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d - cardi\tb:%s c:%s s:%s i:%d p:%d k:%d r:%.2f\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    index += 1
    run_command =  \
    'arguments  = -type cardinality -d AA -a 20 -b %s -c %s -s %s -x 1 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
    output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, r, k, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d - cardi\tb:%s c:%s s:%s i:%d p:%d k:%d r:%.2f x:1\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    index += 1
    run_command =  \
    'arguments  = -type cardinality -d AA -a 20 -b %s -c %s -s %s -x 2 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
    output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, r, k, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d - cardi\tb:%s c:%s s:%s i:%d p:%d k:%d r:%.2f x:2\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    comb2 = [(ci,cp) for ci in cutiter for cp in cutperiter]
    for (ci,cp) in comb2:
        # run command
        index += 1
        run_command =  \
            'arguments  = -type cardinality -d AA -a 20 -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, r, k, index)
        # write to file
        cfile.write(run_command)
        prob_index = '%04d - cardi\tb:%s c:%s s:%s i:%d p:%d k:%d r:%.2f\n' %(index, br, cu, sea, ci, cp, k, r)
        ifile.write(prob_index)



# error      = ./error.txt\n\
# log        = ./log.txt\n\

