CXX = g++ -std=c++1z
CXXFLAGS = -Wall -g 

main: main.o
	$(CXX) $(CSSFLAGS) main.o -o main 

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

clean: 
	rm -f main 
	rm -f *.o 
	rm -f *.dat