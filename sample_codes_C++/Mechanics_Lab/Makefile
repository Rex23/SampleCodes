CXX = g++

INC_DIR = ./include

SOURCE_DIR = ./source

OBJ_DIR = ./source

LAB_DIR = ./lab

CXXFLAGS += -g -Wall -Wextra -pthread -std=c++11 -fpermissive -I$(INC_DIR) -w

SOURCES = $(wildcard $(SOURCE_DIR)/*.cpp)

OBJ = $(SOURCES:.cpp=.o)

MECHANICS_LAB_SRCS_ = $(SOURCES) 

$(LAB_DIR)/LAB: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $^ 

all: $(LAB) 

clean:
	rm -f $(OBJ_DIR)/*.o
