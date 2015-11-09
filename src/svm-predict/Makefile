objs = svm.o svm-predict-ha.o

svm-predict : $(objs)
	g++ -o svm-predict $(objs)

svm.o : svm.cpp svm.h
	g++ -O2 -c svm.cpp

svm-predict-ha.o : svm-predict-ha.cpp
	g++ -O2 -c svm-predict-ha.cpp

clean :
	rm -f svm-predict $(objs)
