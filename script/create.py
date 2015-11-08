
# This create.py script search the data folder and
# create condor submission file (condor.sub) for all data

import os

# Step -1: Remove all files in test portfolio-output folder

directory = '../test/portfolio-output/'
files = os.listdir(directory)
for efile in files:
    if ".txt" in efile:
        os.remove(directory + efile)



# Step 0: Condor for AA problems

# Step 0.0: Arrays
dataset = 'RD1' # 'RD1' or 'AA'
capital = [100000]
ret = [0.06] 
asize = 20 # not yet implemented
cardinality = [5, 10]
branch = ['mf', 'hc'] # 'mf', 'hc', 'bonami', 'hvar', 'random'
cut = ['mf', 'hc', 'bonami', 'hvar'] # 'mf', 'hc', 'bonami', 'hvar', 'random'
search = ['df1'] # 'df0', 'df1', 'best', 'breadth'
cutiter = [1, 6]
cutperiter = [1]
cutlim = [1, 2]
mincutdepth = [0, 3]

# Step 0.1: Open file
cfile = open('condor.sub','w')
ifile = open('exp.txt','w')

# Step 0.2: Print common values
common_command = \
'Executable = ../test/portfolio \n\
Universe   = vanilla\n\
getenv     = true\n\
request_cpus = 1\n\
transfer_executable = false \n\
\n\n'
cfile.write(common_command)

# Step 0.3: Loop over arrays

comb = [(c,r,br,cu,sea) for c in capital for r in ret for br in branch for cu in cut for sea in search]
index = 0

for (c,r,br,cu,sea) in comb:
    index = index + 1
    ci = 0
    cp = 0
    cl = 0
    # run command
    run_command =  \
'arguments  = -type roundlot -d %s -a %s -b %s -c %s -s %s -x 0 -l %d -i %d -p %d -C %d -r %f -f 0 -o 2\n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, cl, ci, cp, c, r, index, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tround\tb:%s\tc:%s\ts:%s\tl:%d\ti:%d\tp:%d\tC:%d\tr:%.4f\tleader\n' %(index, br, cu, sea, cl, ci, cp, c, r)
    ifile.write(prob_index)

    # x = 1
    index += 1
    cl = 20
    run_command =  \
'arguments  = -type roundlot -d %s -a %s -b %s -c %s -s %s -x 1 -l %d -i %d -p %d -C %d -r %f -f 0 -o 2\n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, cl, ci, cp, c, r, index, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tround\tb:%s\tc:%s\ts:%s\tl:%d\ti:%d\tp:%d\tC:%d\tr:%.4f\tx:1\n' %(index, br, cu, sea, cl, ci, cp, c, r)
    ifile.write(prob_index)
    
    comb2 = [(ci,cp,cl,cd) for ci in cutiter for cp in cutperiter for cl in cutlim for cd in mincutdepth]
    for (ci,cp,cl,cd) in comb2:
        # run command
        index += 1
        run_command =  \
'arguments  = -type roundlot -d %s -a %s -b %s -c %s -s %s -x 0 -l %d -i %d -p %d -cd %d -C %d -r %f -f 0 -o 2\n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, cl, ci, cp, cd, c, r, index, index)
        # write to file
        cfile.write(run_command)
        prob_index = '%04d\tround\tb:%s\tc:%s\ts:%s\tl:%d\ti:%d\tp:%d\tC:%d\tr:%.4f\tcd:%d\n' %(index, br, cu, sea, cl, ci, cp, c, r, cd)
        ifile.write(prob_index)

#  -----------  CARDINALITY -----------

comb3 = [(r,br,cu,sea,k) for r in ret for br in branch for cu in cut for sea in search for k in cardinality]
cl = 20
for (r,br,cu,sea,k) in comb3:
    index = index+1
    ci = 0
    cp = 0
    run_command =  \
    'arguments  = -type cardinality -d %s -a %s -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tcardi\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tleader\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    index += 1
    run_command =  \
    'arguments  = -type cardinality -d %s -a %s -b %s -c %s -s %s -x 1 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index,index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tcardi\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tx:1\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    index += 1
    run_command =  \
    'arguments  = -type cardinality -d %s -a %s -b %s -c %s -s %s -x 2 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
    output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tcardi\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tx:2\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    comb2 = [(ci,cp,cd) for ci in cutiter for cp in cutperiter for cd in mincutdepth]
    for (ci,cp,cd) in comb2:
        # run command
        index += 1
        run_command =  \
            'arguments  = -type cardinality -d %s -a %s -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -cd %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, cd, r, k, index, index)
        # write to file
        cfile.write(run_command)
        prob_index = '%04d\tcardi\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tcd:%d\n' %(index, br, cu, sea, ci, cp, k, r, cd)
        ifile.write(prob_index)


#  -----------  SINGLE CARDINALITY -----------

#comb3 = [(r,br,cu,sea,k) for r in ret for br in branch for cu in cut for sea in search for k in cardinality]
for (r,br,cu,sea,k) in comb3:
    index = index+1
    ci = 0
    cp = 0
    run_command =  \
    'arguments  = -type single -d %s -a %s -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tsingle\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tleader\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    index += 1
    run_command =  \
    'arguments  = -type single -d %s -a %s -b %s -c %s -s %s -x 1 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index,index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tsingle\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tx:1\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    index += 1
    run_command =  \
    'arguments  = -type single -d %s -a %s -b %s -c %s -s %s -x 2 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
    output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tsingle\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tx:2\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    index += 1
    ci = 0
    cp = 0
    run_command =  \
    'arguments  = -type single -d %s -a %s, -b %s -c %s -s %s -x 0 -l 0 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct linear \n\
    output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index, index)
    # write to file
    cfile.write(run_command)
    prob_index = '%04d\tsingle\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\tlinear\n' %(index, br, cu, sea, ci, cp, k, r)
    ifile.write(prob_index)
    
    
    comb2 = [(ci,cp) for ci in cutiter for cp in cutperiter]
    for (ci,cp) in comb2:
        # run command
        index += 1
        run_command =  \
            'arguments  = -type single -d %s -a %s -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -r %f -f 0 -o 2 -k %d -ct quadratic \n\
output     = ../test/portfolio-output/out.%04d.txt\
\nlog      = ../test/portfolio-output/log/%04d.txt\
\nqueue 1\n\n' %(dataset, asize, br, cu, sea, ci, cp, r, k, index, index)
        # write to file
        cfile.write(run_command)
        prob_index = '%04d\tsingle\tb:%s\tc:%s\ts:%s\tl:100\ti:%d\tp:%d\tk:%d\tr:%.4f\n' %(index, br, cu, sea, ci, cp, k, r)
        ifile.write(prob_index)






# error      = ./error.txt\n\
# log        = ./log.txt\n\

