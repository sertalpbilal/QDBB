
# This create.py script search the data folder and
# create condor submission file (condor.sub) for all data


# Step 0: Condor for AA problemsOP

# Step 0.0: Arrays
capital = [1000]
ret = [0.04, 0.06]
card = [3, 5, 10]
branch = ['mf','hc','random']
cut = ['mf', 'hc', 'random']
search = ['df0', 'df1', 'breadth', 'best']
cutiter = [0, 1]
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

comb = [(c,r,br,cu,sea,ci,cp) for c in capital for r in ret for br in branch for cu in cut for sea in search for ci in cutiter for cp in cutperiter]
index = 0

for (c,r,br,cu,sea,ci,cp) in comb:
    index = index + 1
    
    # run command
    run_command =  \
'arguments  = -type roundlot -d AA -a 20 -b %s -c %s -s %s -x 0 -l 100 -i %d -p %d -C %d -r %f -f 0 -o 2\n\
output     = ../test/portfolio-output/out.%04d.txt\nqueue 1\n\n' %(br, cu, sea, ci, cp, c, r, index)
    
    # write to file
    cfile.write(run_command)
    prob_index = '%04d - b:%s c:%s s:%s i:%d p:%d C:%d r:%.2f\n' %(index, br, cu, sea, ci, cp, c, r)
    ifile.write(prob_index)


# error      = ./error.txt\n\
# log        = ./log.txt\n\

