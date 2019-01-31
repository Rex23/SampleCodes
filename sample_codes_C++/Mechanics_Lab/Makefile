CXX = g++

INC_DIR = ./include

SOURCE_DIR = ./source

OBJ_DIR = ./source

LAB_DIR = ./lab

BOOST_ROOT = /usr/local/include/boost/

BOOSTFLAGS = -I$(BOOST_ROOT)

CXXFLAGS += -g -Wall -Wextra -pthread -std=c++1y -fpermissive -I$(INC_DIR) -w $(BOOSTFLAGS)

SOURCES = $(wildcard $(SOURCE_DIR)/*.cpp)

OBJ = $(SOURCES:.cpp=.o)

$(LAB_DIR)/LAB: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $^ 

all: $(LAB_DIR)/LAB

#For google test:

GTEST_DIR_LOC = ./googletest

GTEST_DIR = ./googletest

GTEST_LIB_DIR = $(GTEST_DIR)/make

GTEST_LIBS = libgtest.a libgtest_main.a

GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
		$(GTEST_DIR)/include/gtest/internal/*.h

CPPFLAGS += -isystem $(GTEST_DIR)/include
 
$(GTEST_DIR_LOC)/test.o: $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $(GTEST_DIR_LOC)/test.cpp

OBJ_TEST = $(SOURCE_DIR)/Utilities.o
 
$(GTEST_DIR_LOC)/google: $(OBJ_TEST) $(GTEST_DIR_LOC)/test.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L$(GTEST_LIB_DIR) -lgtest_main -lpthread $^ -o $@

google: $(GTEST_DIR_LOC)/google

clean:
	rm -f $(OBJ_DIR)/*.o $(GTEST_DIR_LOC)/*.o

cleangoogle:
	rm -f $(GTEST_DIR_LOC)/*.o
