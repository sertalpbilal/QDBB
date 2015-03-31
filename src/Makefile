
HOME=/usr/local
PACKAGE=/home/sec312/package
SRCDIR=../src
INCLUDES=-I$(HOME)/mosek/7/tools/platform/linux64x86/h \
         -I$(PACKAGE)/lapack-3.5.0/lapacke/include

LIBS=-L$(HOME)/mosek/7/tools/platform/linux64x86/bin \
  -lmosek64 \
  -lblas -llapack

CFLAGS=-ggdb -Wall -Wno-unused-value -g -O0  -Wl,-rpath-link=$(HOME)/mosek/7/tools/platform/linux64x86/bin


DEBUG=0
CC=g++ $(CFLAGS)
SRCFILES=qdbb.cpp \
        ../test/portfolio.cpp
#test_main.cpp
#$(SRCDIR)/portfolio.cpp

  
all: $(SRCFILES)
	$(CC) $(SRCFILES) -o qdbb $(INCLUDES) $(LIBS)
	
#$(CC) -c $(SRCFILES) -o simopt.o
#$(CC) simopt.o busproblem.c -o simopt $(LIBS)
#call problem as follows:
#simopt rspline rsp.opt
#simopt paris option.prs

clean:
	rm -rf *.o *~ qdbb
	#rm -rf $(SRCDIR)/*.o $(SRCDIR)/*~ 




	
