#Tianyang Li, V00814119

CXX=g++
CXXFLAGS=-g -O3 -Wall
LDFLAGS=-lglut -lX11 -lGL -lGLU

all: ECGDemo

ECGDemo: main.o ECGDetect.o Filter.o FatalDetect.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
main.o: main.cpp ECGDetect.hpp FatalDetect.hpp
	$(CXX) -c $(CXXFLAGS) $^
ECGDetect.o: ECGDetect.cpp ECGDetect.hpp Filter.hpp
	$(CXX) -c $(CXXFLAGS) $^
Filter.o: Filter.cpp Filter.hpp
	$(CXX) -c $(CXXFLAGS) $^
FatalDetect.o: FatalDetect.cpp FatalDetect.hpp Filter.hpp
	$(CXX) -c $(CXXFLAGS) $^

clean:
	rm -f *.o ECGDemo *.gch

.PHONY: all clean

