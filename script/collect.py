
# This file gathers the results from test folder and put them into table

import os, os.path, sys, pdb
import string

result = []

status_format  = 0
summary_format = 0
no_cmd = 0
if(len(sys.argv) >= 2):
    if("summary" in sys.argv[1]):
        summary_format = 1
    if("status" in sys.argv[1]):
        status_format = 1
    if("no-cmd" in sys.argv[1]):
        no_cmd = 1

lead_node = 0
lead_time = 0
# outer loop: for all files in system
#directory = '../test/portfolio-output/'
directory = 'output/'
files = os.listdir(directory)
expfile = open('exp.txt','r').readlines()

totalsolved = 0
for efile in sorted(files):
    if ".o" not in efile: #".txt" not in efile:
        continue
    ofile = open(directory + efile, 'r')
    found = 0
    exp = []
    node_id = efile.split('.')[1].lstrip("0")
    #exp.append(node_id)
    for line in ofile:
        words = line.split()
        if "./portfolio" in line and no_cmd == 0:
            #runcmd = line
            runcmd = string.replace(line,'\n','')
            runcmd = string.replace(runcmd,'\r','')
            runcmd = string.replace(runcmd,' ','~')
            exp.append(runcmd)
        if "PID" in line:
            exp.append(words[1])
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
        if "cuts generated:" in line:
            exp.append(words[4])
        if "Total objective improvement:" in line:
            exp.append(words[4])
        if "Best objective improvement:" in line:
            exp.append(words[4])
        if "SOCO solved:" in line:
            exp.append(words[4])
        if "SUMMARY" in line:
            found = 1
            totalsolved += 1
            # exp.append(efile)
    if found:
        '''
        if "leader" in expfile[int(node_id)]:
            exp.append("leader")
            lead_node = exp[2]
            lead_time = exp[5]
            exp.append(" ")
        else:
            node_imp = "0"
            node_imp = "% 8.2f %%" % ( (float(lead_node)-float(exp[2]))/max(1,float(lead_node)) * 100 )
            exp.append(node_imp)
            time_imp = "% 8.2f %%" % ( (float(lead_time)-float(exp[5]))/max(1,float(lead_time)) * 100 )
            exp.append(time_imp)
        time_per_soco = "%.6f" % (float(exp[5])/float(exp[4]))
        exp.append(time_per_soco)
        '''
        disabled_for_now = 1
    if( summary_format == 0 or (summary_format == 1 and found ==1 )):
            result.append(exp)
    ofile.close()

#titles = []
#title = ['ID', 'Filename','\tNProc','NGen','SOCO S','Time       ','Cuts G','Cuts Ap','Total Impr.','Best Impr.', 'Objective','Node Chg','Time Chg','Time/SOCO']
#result.insert(0, title)

print('ID Dataset Asset P.Type CutMethod Branching Cutting Searching Capital Return Risk Cardinality Limit Iter CutPerIter ObjType Filename ID2 Cmd NProc NGen SOCOsolved Time CutsGen CutsApp TotalImp BestImp Objective')

i = 0
if(status_format != 1):
    for ex in result:
        #print ' '.join(ex)
        mex = ' '.join(ex)
        #print(mex)
        try:
            #print('Ex: {}'.format(ex[0]))
            #print('Exp: {}'.format(expfile[i].split(' ', 1)[0]))
            #while(not ex[0] in expfile[i].split(' ', 1)[0]):
            #    i = i+1
            print expfile[i][:-1],
            print mex
            i = i+1
        except IndexError:
            print('IndexError')
            continue

for i in result:
    #print(i)
    a=2
    
# print str(totalsolved) + ' out of ' + str(len(files)-1) + ' experiments are completed.'

#for a in range(1,20):
#    print '%02d' %(a)
