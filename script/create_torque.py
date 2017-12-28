import os
import itertools

problemtype = 2 # 1: roundlot, 2:cardinality, 3:single-bound
file_append = False #True
PID = 1
cps = []

for a in [10, 20,] : #25, 50]: #, 200, 300, 400]:
    for s in range(1,11): #,2,3,4,5,6,7,8,9,10]:
        for N in [50000, 100000]:
            for ret in [0.02, 0.03, 0.04, 0.05, 0.06]:
                for ot in [0]: #[0,1]
                    cps.append(['N_{}_{}'.format(500,s), a, N, ret, ot])
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

search  = ['df1'] #['best', 'df0', 'df1', 'breadth', 'dfc'] # 'best','df0','df1','breadth'] #'best','breadth','df0','df1','dfc']
if problemtype==1:
    cutmeth = [0, 1, 2, 4] # 0: BB, 1: BCC-F, 2: BCC-I, 3: BCC-R, 4: BCC-D, 10: MOSEK, 11: MOSEK-R
elif problemtype==2 or problemtype==3:
    cutmeth = [0, 1, 2, 3]

branch_rule = ['mf']
cut_rule = ['hc']
#branch_rule = ['mf', 'hc', 'bonami', 'hvar'] #['mf'] #,'hc','bonami','hvar'] # 'random'
#cut_rule    = ['mf', 'hc', 'bonami', 'hvar'] # ['mf','hc','bonami','hvar'] # 'random'

if file_append:
    jobs = open('jobs.sh','a')
    jlist = open('exp.txt','a')
else:
    jobs = open('jobs.sh','w+')
    jlist = open('exp.txt','w+')

for b in branch_rule:
    for c in cut_rule:
        for (d,a,C,r,ot) in cps:
            for s in search:
                for x in cutmeth:
                    l = 100
                    # problem specific assignments
                    if problemtype==1:
                        t = 'roundlot'
                    elif problemtype==2:
                        t = 'cardinality'
                        k = round((C+50000)/50000)
                        l = a
                    elif problemtype==3:
                        t = 'single'
                    # BCC-I, 3 iter 1 cut per iter
                    if(x == 2):
                        i = 3
                        p = 1
                    else:
                        i = 1
                        p = 3
                    rsk = r*250
                    allvars = (PID,d,a,t,x,b,c,s,C,r,rsk,k,l,i,p,ot,PID,d,a,r,x)
                    jobs.write('qsub -V -l nodes=1:ppn=1,mem=7GB,vmem=7GB -q batch -v PID=%02d,D=%s,A=%s,T=%s,X=%s,B=%s,C=%s,S=%s,CP=%s,R=%s,RSK=%s,k=%s,L=%s,I=%s,P=%s,OT=%s test.pbs -o output/ -e output/ -N p_%04d_%s_%s_%s_%s\n' % allvars)
                    jlist.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s (p_%04d_%s_%s_%s_%s)\n' % allvars)
                    PID = PID + 1

jobs.close()
jlist.close()

print('Done!')

