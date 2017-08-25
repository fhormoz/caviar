CC=g++
DIC=$(PWD)
CFLAGS=-c -Wall -g  -I $(DIC) 
LDFLAGS= -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl
SOURCES1=caviar.cpp PostCal.cpp Util.cpp TopKSNP.cpp 
SOURCES2=ecaviar.cpp PostCal.cpp Util.cpp
SOURCES3=setcaviar.cpp PostCal.cpp Util.cpp
SOURCES4=multipop.cpp PostCalMultiPop.cpp Util.cpp 
EXECUTABLE1=CAVIAR
EXECUTABLE2=eCAVIAR
EXECUTABLE3=setCAVIAR
EXECUTABLE4=mupCAVIAR

all: $(SOURCES1) $(EXECUTABLE1) $(SOURCES2) $(EXECUTABLE2) $(SOURCES3) $(EXECUTABLE3) $(SOURCES4) $(EXECUTABLE4)
	
$(EXECUTABLE1): $(SOURCES1) 
	$(CC) $(SOURCES1)   $(LDFLAGS) -o $@

$(EXECUTABLE2): $(SOURCES2)
	$(CC) $(SOURCES2)   $(LDFLAGS) -o $@

$(EXECUTABLE3): $(SOURCES3)
	$(CC) $(SOURCES3)   $(LDFLAGS) -o $@

$(EXECUTABLE4): $(SOURCES4)
	$(CC) $(SOURCES4)   $(LDFLAGS) -o $@

clean:
	rm CAVIAR
	rm eCAVIAR
	rm mupCAVIAR
	rm setCAVIAR
