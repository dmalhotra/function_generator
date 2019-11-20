CC = gcc
CXX = g++
F90 = gfortran

CFLAGS+=-I$(CONDA_PREFIX)/include -flto -Wall -g -O3 -ffast-math -march=native
LDFLAGS+=-L$(CONDA_PREFIX)/lib -flto -O3

function_generator.o : function_generator.hpp fg_interface.h fg_interface.cpp
	$(CXX) $(CFLAGS) -std=c++11 -c fg_interface.cpp -o function_generator.o

cpp_example : cpp_example.cpp
	$(CXX) $(CFLAGS) -std=c++11 -c cpp_example.cpp
	$(CXX) cpp_example.o -o cpp_example $(LDFLAGS) -lgsl

c_example : function_generator.o c_example.c
	$(CC) $(CFLAGS) -c c_example.c
	$(CC) c_example.o function_generator.o -lm -lstdc++ -o c_example $(LDFLAGS) -lgsl

fortran_example : function_generator.o fortran_example.f90
	$(F90) $(CFLAGS) -c fg_interface.f90
	$(F90) $(CFLAGS) -c fortran_example.f90
	$(F90) -o fortran_example function_generator.o fg_interface.o fortran_example.o -lstdc++ $(LDFLAGS)

all : cpp_example c_example fortran_example

clean :
	rm -f *.o *.mod cpp_example c_example fortran_example
