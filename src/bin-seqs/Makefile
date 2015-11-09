bin-seqs : bin-seqs.o OftenUsedFun.o seqop.o
	g++ -o bin-seqs bin-seqs.o OftenUsedFun.o seqop.o

bin-seqs.o : bin-seqs.cpp
	g++ -O2 -c bin-seqs.cpp

OftenUsedFun.o : OftenUsedFun.cpp OftenUsedFun.h
	g++ -O2 -c OftenUsedFun.cpp

seqop.o : seqop.cpp seqop.h
	g++ -O2 -c seqop.cpp

clean :
	rm -f bin-seqs bin-seqs.o OftenUsedFun.o seqop.o
