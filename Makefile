CPP    = g++
CFLAGS = -std=c++11 -Wall -O3
LIBFLAGS = -lgsl -lgslcblas -lm -ltinyxml -L/usr/local/lib
INCLUDEFLAGS =-iquote.

SOURCEDIR=./src

# sources
SOURCES = $(shell find $(SOURCEDIR) -name '*.cpp')
OBJECTS = $(SOURCES:%.cpp=%.o)

EXECUTABLE=./bin/inclined

# the whole project
all: $(EXECUTABLE) ./bin/openmp_test ./bin/vertical_test

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(CFLAGS) main.cpp $(OBJECTS) $(INCLUDEFLAGS) -o $@ $(LIBFLAGS)

./bin/openmp_test : ./test/perf/openmp_test.o
	$(CPP) $(CFLAGS) ./test/perf/openmp_test.o $(OBJECTS) $(INCLUDEFLAGS) -o $@ $(LIBFLAGS)

./bin/vertical_test : ./test/vertical_test.o
	$(CPP) $(CFLAGS) ./test/vertical_test.o $(OBJECTS) $(INCLUDEFLAGS) -o $@ $(LIBFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) $(INCLUDEFLAGS) -c $< -o $@ $(LIBFLAGS)

clean:
	rm -f $(OBJECTS) ./test/perf/openmp_test.o ./test/vertical_test.o
	rm -f $(EXECUTABLE) ./bin/openmp_test ./bin/vertical_test
	rm -f *.vts
