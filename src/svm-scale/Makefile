objs = svm.o svm-scale.o

svm-scale : $(objs)
	g++ -o svm-scale $(objs)

svm.o : svm.cpp svm.h
	g++ -O2 -c svm.cpp

svm-scale.o : svm-scale.c
	g++ -O2 -c svm-scale.c

clean :
	rm -f svm-scale $(objs)
