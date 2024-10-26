CXX=g++
CLUSTERC=/DATA/APPS/gcc/7.2.0/bin/g++
COMMON_FLAGS=-Wall -fopenmp
CFLAGS=-O3
FASTFLAGS =-Ofast
OFNAME=currentne2
LIBFILES=lib/*.cpp

currentne:
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -mcmodel=medium -o $(OFNAME) currentne2.cpp $(LIBFILES)
cluster:
	$(CLUSTERC) $(COMMON_FLAGS) $(CFLAGS) -mcmodel=medium -static -o $(OFNAME) currentne2.cpp $(LIBFILES)
static:
	$(CXX) $(COMMON_FLAGS) $(CFLAGS) -mcmodel=medium -static -o $(OFNAME) currentne2.cpp $(LIBFILES)
