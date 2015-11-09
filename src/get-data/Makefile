objs = feature.o get-data.o meta.o OftenUsedFun.o seqop.o tis.o

get-data : $(objs)
	g++ -o get-data $(objs)

feature.o : feature.cpp feature.h
	g++ -O2 -c feature.cpp

get-data.o : get-data.cpp
	g++ -O2 -c get-data.cpp

meta.o : meta.cpp meta.h
	g++ -O2 -c meta.cpp

OftenUsedFun.o : OftenUsedFun.cpp OftenUsedFun.h
	g++ -O2 -c OftenUsedFun.cpp

seqop.o : seqop.cpp seqop.h
	g++ -O2 -c seqop.cpp

tis.o : tis.cpp tis.h
	g++ -O2 -c tis.cpp

clean:
	rm -f get-data $(objs)
