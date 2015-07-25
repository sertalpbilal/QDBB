
HOME=/usr/local
PACKAGE=/home/sec312/package
SRCDIR=../src
INCLUDES=-I$(PACKAGE)/mosek/7/tools/platform/linux64x86/h \
         -I$(PACKAGE)/lapack-3.5.0/lapacke/include

LIBS=-L$(PACKAGE)/mosek/7/tools/platform/linux64x86/bin \
  -lmosek64 \
  -lblas -llapack

CFLAGS=-ggdb -Wall -Wno-unused-value -g -O0  -Wl,-rpath-link=$(PACKAGE)/mosek/7/tools/platform/linux64x86/bin


DEBUG=0
CC=g++ $(CFLAGS)
SRCFILES=qdbb.cpp \
        ../test/portfolio.cpp
SRCCARDINALITY=qdbb.cpp \
        ../test/portfolio.cpp \
        ../test/cardinality.cpp
SRCROUNDLOT=qdbb.cpp \
        ../test/portfolio.cpp \
        ../test/roundlot.cpp
#test_main.cpp
#$(SRCDIR)/portfolio.cpp

all: $(SRCFILES)
	$(CC) $(SRCFILES) -o qdbb $(INCLUDES) $(LIBS)
	# ./qdbb asset rule priority selection limit cut/iter iter C     R    tolerance**
	#              0-5  0-1      0-1       >0    -1, >0   >0              >0
	# ./qdbb 10    1    0        0         100   2        1    10000 0.06 

roundlot: $(SRCROUNDLOT)
	$(CC) $(SRCROUNDLOT) -o roundlot $(INCLUDES) $(LIBS)
            # ./roundlot asset rule priority selection limit cut/iter iter C     R    tolerance**
	#                  0-5  0-1      0-1       >0    -1, >0   >0              >0
	# ./roundlot 10    1    0        0         100   2        1    10000 0.06 

cardinality: $(SRCCARDINALITY)
	$(CC) $(SRCCARDINALITY) -o roundlot $(INCLUDES) $(LIBS)

# argv  #
# 1	# number of variables (asset)
# 2	# rule	0: no cut
# 	#	1: always cut
#	#	2: fading cut
#	#	3: root-heuristic cut
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
# 10    # tolerance: objective tolerance for fading cuts //// **not yet implemented



clean:
	rm -rf *.o *~ qdbb
	#rm -rf $(SRCDIR)/*.o $(SRCDIR)/*~ 




