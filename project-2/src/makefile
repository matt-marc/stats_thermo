CXX = g++ -std=c++1z
CXXFLAGS = -Wall -g 

main: main.o DVS.o
	$(CXX) $(CSSFLAGS) -o main main.o DVS.o

main.o: main.cpp DVS.h
	$(CXX) $(CXXFLAGS) -c main.cpp

DVS.o: DVS.cpp DVS.h
	$(CXX) $(CXXFLAGS) -c DVS.cpp

clean: 
	rm -f main 
	rm -f *.o 
	rm -f *.dat