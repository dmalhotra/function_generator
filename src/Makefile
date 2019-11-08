CXX = g++

FunctionGeneratorPython : function_generator.cpp function_generator.hpp
	$(CXX) -g -O3 -march=native -ffast-math -Wall -shared -std=c++11 -fPIC -DPYTHON_MODULE `python3 -m pybind11 --includes` \
		function_generator.cpp -o FunctionGenerator`python3-config --extension-suffix` -fvisibility=hidden

main : function_generator.hpp main.cpp
	$(CXX) -Wall -g -O3 -ffast-math -march=native -I$(CONDA_PREFIX)/include main.cpp -o main -lgsl

all :
	make FunctionGeneratorPython
	make main

make clean :
	rm -f *.so main
