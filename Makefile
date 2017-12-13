CPP := g++
OPTION := -O3 -ffast-math -funroll-loops -march=native -fopenmp -mavx

SOURCE := 2d_ising.cpp
SOURCE_M := 2d_ising_metropolis.cpp

FILE := ising

FILE:
	$(CPP) $(SOURCE_M) -o $(FILE) $(OPTION)

clean:
	rm -f $(FILE)
