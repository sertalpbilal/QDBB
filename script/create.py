
# This create.py script search the data folder and
# create condor submission file (condor.sub) for all data

import os
import itertools

# Step -1: Remove all files in test portfolio-output folder

directory = '../test/portfolio-output/'
files = os.listdir(directory)
for efile in files:
    if ".txt" in efile:
        os.remove(directory + efile)



# Step 0: Condor for AA problems

# Step 0.0: Arrays
# problems = ['roundlot', 'cardinality', 'single']
dataset = ['RD0', 'RD1'] #, 'RD1', 'RD2', 'RD3', 'RD4', 'RD5', 'RD6', 'RD7', 'RD8', 'RD9'] # 'RD0' to 'RD9' and 'AA'
capital = [200000]
ret = [0.04] 
problemsize = [20] # not yet implemented
cardinality = []
branch = ['mf', 'hc', 'bonami', 'random'] # 'mf', 'hc', 'bonami', 'hvar', 'random'
cut = ['mf', 'hc', 'bonami', 'hvar'] # 'mf', 'hc', 'bonami', 'hvar', 'random'
search = ['df1'] # 'df0', 'df1', 'best', 'breadth'
cutiter = [1]
cutperiter = [1]
cutlim = [1]
mincutdepth = [0]
verbosity = 1
f_ = 0

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
prob_command = \
'arguments  = -type %s -d %s -a %s -b %s -c %s -s %s -x %d -l %d -i %d -p %d -mind %d -C %d -k %d -r %f -ct %s -f 0 -o 1 -dct 0 -maxd 1 \n\
 output     = ../test/portfolio-output/out.%04d.txt\n\
 log      = ../test/portfolio-output/log/%04d.txt\n\
 queue 1\n\n'

# Order: ID, type, dataset, size, branch, cut, search, limit, stratX, iter, cutper, min depth, capital or cardinality, return, quad-linear, leader
log_command = \
'%04d\t\
type:%s\t\
data:%s\t\
n:%d\t\
b:%s\t\
c:%s\t\
s:%s\t\
x:%d\t\
l:%d\t\
i:%d\t\
p:%d\t\
mind:%d\t\
C:%d\t\
k:%d\t\
r:%.4f\t\
ct:%s\t\
info:%s\n'
ifile.write('Problem Definitions\n')

allproblems = []
comb = list(itertools.product(cutlim,cutiter,cutperiter,mincutdepth))
index = 0;

# Problem 1 - Roundlot
pname = 'roundlot'
for (ds,ps,br,cu,sea,c,r) in list(itertools.product(dataset,problemsize,branch,cut,search,capital,ret)):
    x = cl = ci = cp = cd = k = 0
    index, ct, pinfo = index+1, 'null', 'leader'
    temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity,pinfo)
    allproblems.append(temp)
    x = 2
    for cl in cutlim:
        index, x, pinfo = index+1, 2, ''
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)
    x = 1
    for (cl,ci,cp,cd) in comb:
        index = index+1
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)
    x = 4
    index, ct, pinfo = index+1, 'null', 'depthorder'
    cl, ci, cp = ps, 1, ps
    temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity,pinfo)
    allproblems.append(temp)

# Problem 2 - Cardinality
pname = 'cardinality'
for (ds,ps,br,cu,sea,k,r) in list(itertools.product(dataset,problemsize,branch,cut,search,cardinality,ret)):
    x = cl = ci = cp = cd  = 0
    index, ct, pinfo = index+1, 'quadratic', 'leader'
    temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity,pinfo)
    allproblems.append(temp)
    index, ct, pinfo = index+1, 'linear', ''
    temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity,pinfo)
    allproblems.append(temp)
    for cl in cutlim:
        index, ct, x, pinfo = index+1, 'quadratic', 2, ''
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)
    for cl in cutlim:
        index, x, pinfo = index+1, 3, ''
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)
    x = 1
    for (cl,ci,cp,cd) in comb:
        index = index+1
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)

# Problem 3 - Single Bound
pname = 'single'
for (ds,ps,br,cu,sea,k,r) in list(itertools.product(dataset,problemsize,branch,cut,search,cardinality,ret)):
    x = cl = ci = cp = cd  = 0
    index, ct, pinfo = index+1, 'quadratic', 'leader'
    temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity,pinfo)
    allproblems.append(temp)
    index, ct, pinfo = index+1, 'linear', ''
    temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity,pinfo)
    allproblems.append(temp)
    for cl in cutlim:
        index, ct, x, pinfo = index+1, 'quadratic', 2, ''
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)
    for cl in cutlim:
        index, x, pinfo = index+1, 3, ''
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)
    x = 1
    for (cl,ci,cp,cd) in comb:
        index = index+1
        temp = (index,pname,ds,ps,br,cu,sea,x,cl,ci,cp,cd,c,k,r,ct,f_,verbosity, pinfo)
        allproblems.append(temp)

# Problem 4 - Combined



for p in allproblems:
    run_command = prob_command %(p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[0],p[0] )
    cfile.write(run_command)
    pinfo = 'leader'
    prob_index = log_command %(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[18])
    ifile.write(prob_index)
            


quit()


