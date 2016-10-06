main: main.o nbs.o
	g++ -std=c++14 -fopenmp -O3 -o main main.o nbs.o

main.o: main.cpp
	g++ -std=c++14 -fopenmp -O3 -c main.cpp

nbs.o: nbs.cpp
	g++ -std=c++14 -fopenmp -O3 -c nbs.cpp nbs.h 

clean:
	rm -f nbs.o nbs.h.gch main.o main
