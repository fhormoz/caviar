CC=g++
CFLAGS=-c -Wall -g -O3
LDFLAGS=-I/usr/local/include/ -L/usr/local/lib/ -lgsl -lgslcblas -lm
SOURCES=caviar.cpp PostCal.cpp Util.cpp TopKSNP.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=CAVIAR

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
	rm *.o	

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm CAVIAR
