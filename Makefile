CC = g++ -Wall -w -std=c++0x -Wl,--no-as-needed -Ofast -march=native

LSOURCE =  utils.cpp seq_model.cpp kmer.cpp dist_model.cpp output.cpp main.cpp
LHEADER =  utils.h SimpleMatrix.h seq_model.h kmer.h dist_model.h output.h

UltASeqAn: $(LSOURCE) $(HEADER)
	$(CC)  $(LSOURCE) -pthread -o UltASeqAn

clean:
	rm -f *.o UltASeqAn
