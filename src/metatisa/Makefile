objs = GeneInfo.o GeneSeq.o Mainprocess.o Metagenome.o OftenUsedFun.o \
	   OftenUsedOperatLib.o Operate.o SequenceTransform.o Trans.o

metatisa : $(objs)
	g++ -o metatisa $(objs)

GeneInfo.o : GeneInfo.cpp GeneInfo.h
	g++ -O2 -c GeneInfo.cpp

GeneSeq.o : GeneSeq.cpp GeneSeq.h
	g++ -O2 -c GeneSeq.cpp

Mainprocess.o : Mainprocess.cpp
	g++ -O2 -c Mainprocess.cpp

Metagenome.o : Metagenome.cpp Metagenome.h
	g++ -O2 -c Metagenome.cpp

OftenUsedFun.o : OftenUsedFun.cpp OftenUsedFun.h
	g++ -O2 -c OftenUsedFun.cpp

OftenUsedOperatLib.o : OftenUsedOperatLib.cpp OftenUsedOperatLib.h
	g++ -O2 -c OftenUsedOperatLib.cpp

Operate.o : Operate.cpp Operate.h
	g++ -O2 -c Operate.cpp

SequenceTransform.o : SequenceTransform.cpp SequenceTransform.h
	g++ -O2 -c SequenceTransform.cpp

Trans.o : Trans.cpp Trans.h
	g++ -O2 -c Trans.cpp

clean :
	rm -f metatisa $(objs)
