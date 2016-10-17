
HOME=/usr/local
PACKAGE=/home/sec312/package
SRCDIR=$(PACKAGE)/qdbb/src
TESTDIR=$(PACKAGE)/qdbb/test
INCLUDES=-I$(PACKAGE)/mosek/7/tools/platform/linux64x86/h \
         -I$(PACKAGE)/lapack-3.5.0/lapacke/include

LIBS=-L$(PACKAGE)/mosek/7/tools/platform/linux64x86/bin \
  -lmosek64 \
  -lblas -llapack

CFLAGS=-ggdb -Wall -Wunused-variable -g -O0  -Wl,-rpath-link=$(PACKAGE)/mosek/7/tools/platform/linux64x86/bin


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




