CC = gcc
CXX = g++
F90 = gfortran

libfunction_generator.so : function_generator.hpp fg_interface.h fg_interface.cpp
	$(CXX) -Wall -g -O3 -ffast-math -march=native -std=c++11 -shared -fPIC -I$(CONDA_PREFIX)/include -o libfunction_generator.so fg_interface.cpp -fpermissive

cpp_example : libfunction_generator.so cpp_example.cpp
	$(CXX) -Wall -g -O3 -ffast-math -march=native -std=c++11 -I$(CONDA_PREFIX)/include -c cpp_example.cpp
	$(CXX) cpp_example.o -o cpp_example -L$(CONDA_PREFIX)/lib -Wl,-rpath $(CONDA_PREFIX)/lib -lgsl -Wl,-rpath $(PWD) -L./ -lfunction_generator

c_example : libfunction_generator.so c_example.c
	$(CC) -Wall -g -O3 -ffast-math -march=native -I$(CONDA_PREFIX)/include -c c_example.c
	$(CC) c_example.o -lm -lstdc++ -o c_example -Wl,-rpath $(PWD) -L./ -lfunction_generator -lgsl

fortran_example : libfunction_generator.so fortran_example.f90
	$(F90) -Wall -g -O3 -ffast-math -march=native -c fg_interface.f90
	$(F90) -Wall -g -O3 -ffast-math -march=native -c fortran_example.f90
	$(F90) fg_interface.o fortran_example.o -o fortran_example -Wl,-rpath $(PWD) -L./ -lfunction_generator

all :
	make cpp_example c_example fortran_example

clean :
	rm -f *.so *.o *.mod cpp_example c_example fortran_example
