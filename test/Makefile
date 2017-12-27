
HOME=/usr/local
PACKAGE=/home/sec312/package
SRCDIR=$(PACKAGE)/qdbb/src
TESTDIR=$(PACKAGE)/qdbb/test
MSKVER=8
INCLUDES=-I$(PACKAGE)/mosek/${MSKVER}/tools/platform/linux64x86/h \
         -I$(PACKAGE)/lapack-3.5.0/lapacke/include

LIBS=-L$(PACKAGE)/mosek/${MSKVER}/tools/platform/linux64x86/bin \
  -lmosek64 \
  -lblas -llapack

CFLAGS=-ggdb -Wall -Wunused-variable -g -O0  -Wl,-rpath-link=$(PACKAGE)/mosek/${MSKVER}/tools/platform/linux64x86/bin

DEBUG=0
CC=g++ $(CFLAGS)
SRCFILES=$(SRCDIR)/qdbb.cpp \
         $(TESTDIR)/portfolio.cpp
SRCCARDINALITY=$(SRCDIR)/qdbb.cpp \
        $(TESTDIR)/portfolio.cpp \
        $(TESTDIR)/cardinality.cpp
SRCROUNDLOT=$(SRCDIR)/qdbb.cpp \
        $(TESTDIR)/portfolio.cpp \
        $(TESTDIR)/roundlot.cpp
SRCPORT = $(SRCDIR)/qdbb.cpp \
	$(TESTDIR)/portfolio.cpp \
	$(TESTDIR)/roundlot.cpp \
	$(TESTDIR)/cardinality.cpp \
	$(TESTDIR)/single.cpp \
	$(TESTDIR)/diverse.cpp \
	$(TESTDIR)/combined.cpp
#test_main.cpp
#$(SRCDIR)/portfolio.cpp

all: $(SRCFILES)
	$(CC) $(SRCFILES) -o qdbb $(INCLUDES) $(LIBS)
	# ./qdbb asset rule priority selection limit cut/iter iter C     R    tolerance**
	#              0-5  0-1      0-1       >0    -1, >0   >0              >0
	# ./qdbb 10    1    0        0         100   2        1    10000 0.06 

portfolio: $(SRCPORT)
	$(CC) $(SRCPORT) -o portfolio $(INCLUDES) $(LIBS)
	# ./portfolio -type          -d    -a    -C        -r      -b    -c    -s     -l     -x   -i   -p   -o   -f
	# ./portfolio -type roundlot -d AA -a 10 -C 100000 -r 0.06 -b mf -c mf -s df0 -l 100 -x 0 -i 1 -p 2 -o 1 -f 0
	# ./portfolio -t cardinality -d AA -a  3 -k 2      -r 0.01 -b mf -c hc -s df1 -l 100 -x 0 -i 1 -p 2 -o 3 -f 0

roundlot: $(SRCROUNDLOT)
	$(CC) $(SRCROUNDLOT) -o roundlot $(INCLUDES) $(LIBS)
	# ./roundlot asset rule priority selection limit cut/iter iter C     R    tolerance**
	#                  0-5  0-1      0-1       >0    -1, >0   >0              >0
	# ./roundlot 10    1    0        0         100   2        1    10000 0.06 

cardinality: $(SRCCARDINALITY)
	$(CC) $(SRCCARDINALITY) -o cardinal $(INCLUDES) $(LIBS)
	# ./cardinal asset rule priority selection limit cut/iter iter C     R     k   cardinalitytype    tolerance**
	#                  0-5  0-1      0-1       >0    -1, >0   >0                   0 or 1 or 2        >0
	# ./cardinal 10    1    0        0         100   2        1    10000 0.06  4   2

# argv  #
# 1	# number of variables (asset)
# 2	# rule	0: no cut
# 	#	1: always cut
#	#	2: fading cut
#	#	3: root-heuristic cut
#           #           4: not implemented
#           #           5: all cuts at root node (exhaustive binary!)
# 3	# priority	0: most fractional
#	#		1: best return value
# 4	# selection	0: always cut
# 	#		1: only-deep cut
# 5	# limit		maximum allowed cut (for rule 1 and 2)
# 6	# cut/iter	cut per iteration (for rule 1)
#	#		-1: no limit
# 7	# iter		maximum cut generation iteration
# 8	# C		capital
# 9	# R		expected return
# 10    # k             cardinality number
# 11    # tolerance: objective tolerance for fading cuts //// **not yet implemented
# 11    # cardinality type: 0: QC, 1: LC, 2: CONIC


clean:
	rm -rf *.o *~ qdbb
	rm -rf ../test/*.int ../test/*.mps ../test/*.lp 
	rm -rf result/*
	rm portfolio
	#rm -rf $(SRCDIR)/*.o $(SRCDIR)/*~ 


# options for 'portfolio'
# Format:
# -option -alternative info (values,default*)
#
# Mandatory:
#  -type -t Problem type     (roundlot, cardinality, single, diverse, combined)
#  -data -d Data set         (AA,RD0-RD9)
#  -v    -a Number of assets
#  -C       Capital
#  -r       Return rate
#  -k       Cardinality
#  -qa      Cardinality(quad)
#  -ct      Card. type       (quadratic*,linear)
#  
# Optional:
#  -o           Output verbosity (1-6,1*)
#  -f           File output      (0*,1)
#  -b           Branching rule   (mf*,hc,random,bonami,hvar)
#  -c           Cut rule         (mf,hc*,random,bonami,hvar)
#  -s           Search rule      (df0,df1*,breadth,best,dfc)
#  -l           Cut limit        (0-1000,100*)
#  -i           Cut gen. iter    (0-100,1*)
#  -p           Cut per iter     (0-100,3*)
#  -x           Cut strategy     (0 [BB], 1 [BCC-F], 2* [BCC-I], 3 [BCC-R], 4 [BCC-D], 10 [MOSEK], 11 [MOSEK-R])
#  -mind        Min cut depth    (0-100,0*)
#  -maxd        Max cut depth    (0-100,100*)
#  -dct         Deep cut thresh  (0-1,0*)
#  -timelimit   Time limit       (0-1000,1000*)
#  -dryrun  -dr Dry run          (0*,1)
#  -objtype -ot Objective Type  (0*,1,2) min risk, max return, risk-averse



