objs = metalocs-operate.o OftenUsedFun.o

metalocs-operate : $(objs)
	g++ -o metalocs-operate $(objs)

metalocs-operate.o : metalocs-operate.cpp
	g++ -O2 -c metalocs-operate.cpp

OftenUsedFun.o : OftenUsedFun.cpp OftenUsedFun.h
	g++ -O2 -c OftenUsedFun.cpp

clean :
	rm -f metalocs-operate $(objs)
