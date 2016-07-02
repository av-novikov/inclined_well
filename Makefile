CPP    = mpic++
CFLAGS = --std=c++11 -Wall -g
LIBFLAGS = -fopenmp -lgsl -lgslcblas -lm -ltinyxml -lvtkCommonCore-7.0 -lvtkFiltersCore-7.0 -lvtkIOCore-7.0 -lvtkIOXML-7.0 -lvtkCommonDataModel-7.0 -L/usr/local/lib
INCLUDEFLAGS =-iquote. -I/usr/local/include/vtk-7.0

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
	rm -f $(OBJECTS) ./test/perf/2dVS3d_test.o ./test/vertical_test.o
	rm -f $(EXECUTABLE) ./bin/2dVS3d_test ./bin/vertical_test
	rm -f *.vts
