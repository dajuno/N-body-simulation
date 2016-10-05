main: main.o NBS.o
	g++ -Wall -std=c++14 -fopenmp -O3 -o main main.o NBS.o

main.o: main.cpp
	g++ -Wall -std=c++14 -fopenmp -O3 -c main.cpp

NBS.o: NBS.cpp
	g++ -Wall -std=c++14 -fopenmp -O3 -c NBS.cpp NBS.h 

clean:
	rm -f NBS.o main.o main
