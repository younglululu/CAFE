CC = g++ -Wall -w -std=c++0x -Wl, -Ofast -march=native

LSOURCE =  utils.cpp seq_model.cpp kmer.cpp dist_model.cpp output.cpp main.cpp
LHEADER =  utils.h SimpleMatrix.h seq_model.h kmer.h dist_model.h output.h

cafe: $(LSOURCE) $(HEADER)
	$(CC)  $(LSOURCE) -pthread -o cafe

clean:
	rm -f *.o cafe
