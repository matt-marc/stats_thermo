CXX = g++ -std=c++1z
CXXFLAGS = -Wall -g

all: q1 q2

q1: q1.o DVS.o
	$(CXX) $(CSSFLAGS) -o q1 q1.o DVS.o

q1.o: q1.cpp DVS.h
	$(CXX) $(CXXFLAGS) -c q1.cpp

q2: q2.o DVS.o
	$(CXX) $(CSSFLAGS) -o q2 q2.o DVS.o

q2.o: q2.cpp DVS.h
	$(CXX) $(CXXFLAGS) -c q2.cpp

DVS.o: DVS.cpp DVS.h
	$(CXX) $(CXXFLAGS) -c DVS.cpp

main.o: main.cpp DVS.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f main
	rm -f q1
	rm -f q2
	rm -f *.o
	rm -f *.dat
