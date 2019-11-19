CXX = g++

shared : function_generator.hpp function_generator.hpp
	$(CXX) -Wall -g -O3 -ffast-math -march=native -std=c++11 -shared -fPIC -I$(CONDA_PREFIX)/include -o libfunction_generator.so function_generator.cpp fg_interface.cpp -fpermissive

main : function_generator.hpp main.cpp shared
	$(CXX) -Wall -g -O3 -ffast-math -march=native -std=c++11 -I$(CONDA_PREFIX)/include -c main.cpp
	$(CXX) main.o -o main -L$(CONDA_PREFIX)/lib -Wl,-rpath $(CONDA_PREFIX)/lib -lgsl -Wl,-rpath $(PWD) -L./ -lfunction_generator

test : function_generator.hpp test.c fg_interface.cpp function_generator.cpp
	$(CC) -Wall -g -O3 -ffast-math -march=native -I$(CONDA_PREFIX)/include -c test.c
	$(CC) fg_interface.o test.o -lm -lstdc++ -o test -Wl,-rpath $(PWD) -L./ -lfunction_generator -lgsl


all :
	make main

clean :
	rm -f *.so main
