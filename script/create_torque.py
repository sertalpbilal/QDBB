import os
import itertools

file_append = False #True
PID = 1
cps = []

if 0:
    cps =[
        ['N_100_1', 100, 200000, 0.05],
        ['N_100_1', 100, 200000, 0.07],
        ['N_100_1', 100, 400000, 0.09],
        ['N_100_2', 100, 100000, 0.05],
        ['N_100_2', 100, 200000, 0.07],
        ['N_100_2', 100, 200000, 0.09],
        ['N_100_3', 100, 200000, 0.05],
        ['N_100_3', 100, 200000, 0.07],
        ['N_100_3', 100, 200000, 0.09],
        ['N_100_4', 100, 200000, 0.05],
        ['N_100_4', 100, 200000, 0.07],
        ['N_100_4', 100, 200000, 0.09],
        ['N_100_6', 100, 100000, 0.05],
        ['N_100_6', 100, 100000, 0.09],
        ['N_100_6', 100, 200000, 0.05],
        ['N_100_6', 100, 200000, 0.07],
        ['N_100_6', 100, 400000, 0.09],
        ['N_100_7', 100, 100000, 0.05],
        ['N_100_7', 100, 200000, 0.07],
        ['N_100_7', 100, 200000, 0.09],
        ['N_100_8', 100, 100000, 0.05],
        ['N_100_8', 100, 200000, 0.07],
        ['N_100_8', 100, 200000, 0.09],
        ['N_100_9', 100, 200000, 0.05],
        ['N_100_9', 100, 200000, 0.07],
        ['N_100_10', 100, 100000, 0.05],
        ['N_100_10', 100, 100000, 0.07],
        ['N_100_10', 100, 200000, 0.05],
        ['N_100_10', 100, 200000, 0.07],
        ['N_100_10', 100, 200000, 0.09],
        ['N_100_10', 100, 400000, 0.09],
        ['N_100_5', 100, 100000, 0.02]
    ]
else:
    for a in [25, 50, 75, 100]: #, 200, 300, 400]:
        for s in [1,2,3,4,5,6,7,8,9,10]:
            for N in [50000, 100000, 150000, 200000]:
                for ret in [0.02, 0.03, 0.04, 0.05]:
                    cps.append(['N_{}_{}'.format(500,s), a, N, ret])
#cps.append(['AA' ,20,400000, 0.05])
#cps.append(['RD0',10,400000, 0.08])
#cps.append(['RD1',30,100000, 0.07])
#cps.append(['RD2',10,300000, 0.03])
#cps.append(['RD3',10,300000, 0.04])
#cps.append(['RD4',10,800000, 0.03])
#cps.append(['RD5',10,700000,0.05])
#cps.append(['RD6',10,600000, 0.02])
#cps.append(['RD7',10,1000000, 0.05])
#cps.append(['RD8',10,400000, 0.07])
#cps.append(['RD9',10,700000, 0.08])

search  = ['best'] #'best','breadth','df0','df1']
cutmeth = [0, 1, 2, 4] # 0: BB, 1: BCC-F, 2: BCC-I, 3: BCC-R, 4: BCC-D, 10: MOSEK, 11: MOSEK-R

branch_rule = ['mf'] #,'hc','bonami','hvar'] # 'random'
cut_rule    = ['hc'] # ['mf','hc','bonami','hvar'] # 'random'

if file_append:
    jobs = open('jobs.sh','a')
    jlist = open('exp.txt','a')
else:
    jobs = open('jobs.sh','w+')
    jlist = open('exp.txt','w+')

for b in branch_rule:
    for c in cut_rule:
        for (d,a,C,r) in cps:
            for s in search:
                for x in cutmeth:
                    t = 'roundlot'
                    l = 100
                    if(cutmeth == 2):
                        i = 3
                        p = 1
                    else:
                        i = 1
                        p = 3
                    allvars = (PID,d,a,t,x,b,c,s,C,r,l,i,p,PID,d,a,r,x)
                    jobs.write('qsub -V -l nodes=1:ppn=1,mem=7GB,vmem=7GB -q long -v PID=%02d,D=%s,A=%s,T=%s,X=%s,B=%s,C=%s,S=%s,CP=%s,R=%s,L=%s,I=%s,P=%s test.pbs -o output/ -e output/ -N p_%04d_%s_%s_%s_%s\n' % allvars)
                    jlist.write('%s %s %s %s %s %s %s %s %s %s %s %s %s (p_%04d_%s_%s_%s_%s)\n' % allvars)
                    PID = PID + 1

jobs.close()
jlist.close()

print('Done!')

