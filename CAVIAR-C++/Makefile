CC=g++
DIC=$(PWD)
CFLAGS=-c -Wall -g -O3  -I $(DIC) 
LDFLAGS= -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl
SOURCES1=caviar.cpp PostCal.cpp Util.cpp TopKSNP.cpp 
SOURCES2=ecaviar.cpp PostCal.cpp Util.cpp
EXECUTABLE1=CAVIAR
EXECUTABLE2=eCAVIAR

all: $(SOURCES1) $(EXECUTABLE1) $(SOURCES2) $(EXECUTABLE2)
	
$(EXECUTABLE1): $(SOURCES1) 
	$(CC) $(SOURCES1)   $(LDFLAGS) -o $@

$(EXECUTABLE2): $(SOURCES2)
	$(CC) $(SOURCES2)   $(LDFLAGS) -o $@

clean:
	rm CAVIAR
	rm eCAVIAR
