parse-blast : parse-blast.o
	g++ -o parse-blast parse-blast.o

parse-blast.o : parse-blast.cpp
	g++ -O2 -c parse-blast.cpp

clean:
	rm -f parse-blast parse-blast.o
