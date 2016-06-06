CPP    = mpic++
CFLAGS = --std=c++11 -Wall -g
LIBFLAGS = -lgsl -lgslcblas -lm -ltinyxml -lvtkCommonCore-6.0 -lvtkFiltersCore-6.0 -lvtkIOCore-6.0 -lvtkIOXML-6.0 -lvtkCommonDataModel-6.0 -L/usr/local/lib
INCLUDEFLAGS =-iquote. -I/usr/include/vtk-6.0

SOURCEDIR=./src
TESTSDIR=./test

# sources
SOURCES = $(shell find $(SOURCEDIR) -name '*.cpp')
OBJECTS = $(SOURCES:%.cpp=%.o)

# tests
TESTS =$(shell find $(TESTSDIR) -name '*.cpp')
TESTS_OBJ = $(TESTS:%.cpp=%.o)

# the whole project

EXECUTABLE=./bin/inclined
EX_TESTS=./bin/2dVS3d

all: $(EXECUTABLE) $(EX_TESTS)

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(CFLAGS) main.cpp $(OBJECTS) $(INCLUDEFLAGS) -o $@ $(LIBFLAGS)

$(EX_TESTS): $(TEST_OBJ)
	$(CPP) $(CFLAGS) $(TESTS) $(OBJECTS) $(INCLUDEFLAGS) -o $@ $(LIBFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) $(INCLUDEFLAGS) -c $< -o $@ $(LIBFLAGS)

clean:
	rm -f $(OBJECTS) $(TESTS_OBJ)
	rm -f $(EXECUTABLE) $(EX_TESTS)
	rm -f *.vts
