CXX = g++

CXXFLAGS += -O2 -Wall -std=c++17

SOURCES = $(wildcard *.cpp)

OBJ = $(SOURCES:.cpp=.o)

Matching_Engine: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm *.o
