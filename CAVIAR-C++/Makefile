CC=g++
DIC=$(PWD)
CFLAGS=-c -Wall -g -O3  -I $(DIC) 
LDFLAGS= -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl
SOURCES=caviar.cpp PostCal.cpp Util.cpp TopKSNP.cpp 
EXECUTABLE=CAVIAR

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(SOURCES) 
	$(CC) $(SOURCES)   $(LDFLAGS) -o $@

clean:
	rm CAVIAR
