CPP    = mpic++
CFLAGS = --std=c++11 -Wall -g
LIBFLAGS = -lgsl -lgslcblas -lm -ltinyxml -lvtkCommonCore-7.0 -lvtkFiltersCore-7.0 -lvtkIOCore-7.0 -lvtkIOXML-7.0 -lvtkCommonDataModel-7.0 -L/usr/local/lib
INCLUDEFLAGS =-I/usr/local/include/vtk-7.0

SOURCES = $(shell find $(SOURCEDIR) -name '*.cpp')
OBJECTS = $(SOURCES:%.cpp=%.o)

# the whole project

EXECUTABLE=inclined

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(CFLAGS) $(OBJECTS) -o $@ $(LIBFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) $(INCLUDEFLAGS) -c $< -o $@ $(LIBFLAGS)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXECUTABLE)
	rm -f *.vts
