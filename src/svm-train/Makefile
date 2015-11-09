objs = svm.o svm-train-CDS.o
svm-train : $(objs)
	g++ -o svm-train $(objs)

svm.o : svm.cpp svm.h
	g++ -O2 -c svm.cpp

svm-train-CDS.o : svm-train-CDS.cpp
	g++ -O2 -c svm-train-CDS.cpp

clean :
	rm -f svm-train $(objs)
