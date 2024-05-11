CXX = g++
CFLAGS = -fopenmp -O3

main: main.cpp
	$(CXX) -g $(CFLAGS) -o main main.cpp

clean:
	rm -f main
