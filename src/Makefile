
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
#test_main.cpp
#$(SRCDIR)/portfolio.cpp

  
all: $(SRCFILES)
	$(CC) $(SRCFILES) -o qdbb $(INCLUDES) $(LIBS)
        # ./qdbb asset rule priority selection limit cut/iter iter C     R
	#              0-5  0-1      0-1       >0    -1, >0   >0
	# ./qdbb 10    1    0        0         100   2        1    10000 0.06

# rule	0: no cut
#	1: always cut
#	2: fading cut
#	3: root-heuristic cut
# priority	0: most fractional
#		1: best return value
# selection	0: always cut
#		1: only-deep cut
# limit		maximum allowed cut (for rule 1 and 2)
# cut/iter	cut per iteration (for rule 1)
#		-1: no limit
# iter		maximum cut generation iteration
# C		capital
# R		expected return



#$(CC) -c $(SRCFILES) -o simopt.o
#$(CC) simopt.o busproblem.c -o simopt $(LIBS)
#call problem as follows:
#simopt rspline rsp.opt
#simopt paris option.prs

clean:
	rm -rf *.o *~ qdbb
	#rm -rf $(SRCDIR)/*.o $(SRCDIR)/*~ 




