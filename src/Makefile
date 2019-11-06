FunctionGeneratorPython : function_generator.cpp function_generator.hpp
	c++ -g -O3 -march=native -ffast-math -Wall -shared -std=c++11 -fPIC -DPYTHON_MODULE `python3 -m pybind11 --includes` \
		function_generator.cpp -o FunctionGenerator`python3-config --extension-suffix` -fvisibility=hidden

libFunctionGenerator.so : function_generator.cpp function_generator.hpp
	c++ -I$(CONDA_PREFIX)/include -g -O3 -march=native -ffast-math -Wall -shared -std=c++11 -fPIC function_generator.cpp -o libFunctionGenerator.so

main : function_generator.cpp function_generator.hpp main.cpp libFunctionGenerator.so
	make libFunctionGenerator.so
	c++ -O2 -Wall -g -ffast-math -march=native -I$(CONDA_PREFIX)/include main.cpp libFunctionGenerator.so -o main

all :
	make FunctionGeneratorPython
	make libFunctionGenerator.so
	make main

make clean :
	rm -f *.so main
