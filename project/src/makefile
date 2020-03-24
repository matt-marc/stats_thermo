CXX = g++ -std=c++1z
CXXFLAGS = -Wall -g 

main: main.o DVS.o Domain.o
	$(CXX) $(CSSFLAGS) -o main main.o DVS.o Domain.o

main.o: main.cpp DVS.h Domain.h
	$(CXX) $(CXXFLAGS) -c main.cpp

DVS.o: DVS.cpp DVS.h Domain.h
	$(CXX) $(CXXFLAGS) -c DVS.cpp

Domain.o: Domain.cpp Domain.h
	$(CXX) $(CXXFLAGS) -c Domain.cpp

clean: 
	rm -f main 
	rm -f *.o 
	rm -f *.dat