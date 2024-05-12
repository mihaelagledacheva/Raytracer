CXX = g++
CFLAGS = -fopenmp -O3

main: main.cpp
	$(CXX) -g $(CFLAGS) -o main main.cpp

image: image.cpp
	$(CXX) -g $(CFLAGS) -o image image.cpp

clean:
	rm -f main
	rm -f image
