CC=g++
CFLAGS=-O2
LFLAGS=-lfftw3 -lpng

main: main.cpp
	g++ main.cpp -o m $(LFLAGS) $(CFLAGS)