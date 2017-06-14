import os
import itertools

cps = []
cps.append(['AA' ,20,400000, 0.05])
cps.append(['RD0',10,400000, 0.08])
cps.append(['RD1',30,100000, 0.07])
cps.append(['RD2',10,300000, 0.03])
cps.append(['RD3',10,300000, 0.04])
cps.append(['RD4',10,800000, 0.03])
cps.append(['RD5',10,700000,0.05])
cps.append(['RD6',10,600000, 0.02])
cps.append(['RD7',10,1000000, 0.05])
cps.append(['RD8',10,400000, 0.07])
cps.append(['RD9',10,700000, 0.08])

PID = 1
search  = ['best'] #'best','breadth','df0','df1']
cutmeth = [0, 1, 2, 4] # 0: BB, 1: BCC-F, 2: BCC-I, 3: BCC-R, 4: BCC-D, 10: MOSEK, 11: MOSEK-R

branch_rule = ['mf','hc','bonami','hvar'] # 'random'
cut_rule    = ['mf','hc','bonami','hvar'] # 'random'

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
                    allvars = (PID,d,a,t,x,b,c,s,C,r,l,i,p)
                    jobs.write('qsub -V -l nodes=1,mem=4GB,vmem=4GB -q verylong -v PID=%s,D=%s,A=%s,T=%s,X=%s,B=%s,C=%s,S=%s,CP=%s,R=%s,L=%s,I=%s,P=%s test.pbs -o output/ -e output/\n' % allvars)
                    jlist.write('%s %s %s %s %s %s %s %s %s %s %s %s %s\n' % allvars)
                    PID = PID + 1

jobs.close()
jlist.close()

print('Done!')

