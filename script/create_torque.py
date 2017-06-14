import os
import itertools

cps = []
cps.append(['AA' ,20,400000, 0.05])
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

search  = ['best','breadth','df0','df1']
cutmeth = [0, 1, 2, 3] # 0: BB, 1: BCC-F, 2: BCC-I, 3: BCC-R, 4: BCC-D, 10: MOSEK, 11: MOSEK-R

branch_rule = ['mf','hc','random','bonami','hvar']
cut_rule    = ['mf','hc','random','bonami','hvar']

jobs = open('jobs.sh','w+')
jlist = open('exp.txt','w+')

for (d,a,C,r) in cps:
    for s in search:
        for x in cutmeth:
            t = 'roundlot'
            l = 100
            i = 3
            p = 3
            allvars = (d,a,t,x,s,C,r,l,i,p)
            jobs.write('qsub -V -l nodes=1,mem=4GB,vmem=4GB -q verylong -v D=%s,A=%s,T=%s,X=%s,S=%s,C=%s,R=%s,L=%s,I=%s,P=%s test.pbs -o output/ -e output/\n' % allvars)
            jlist.write('%s %s %s %s %s %s %s %s %s %s\n' % allvars) 

jobs.close()
jlist.close()

print('Done!')

