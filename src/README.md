# Experimental C++ variant of function_generator

Build is very low dependency. Should only require a python3 install `Eigen` and `pybind11`.
If not using conda for eigen, you might have to modify `Makefile` to your `Eigen` install path.
```
conda install pybind11 eigen
make all
```
will build a test program `main`, a C++ shared object `libFunctionGenerator.so`, and the python module.

To test the python module (see the `test.py` file included):
```
python3 test.py
```

To test the C++ module (see the `main.cpp` file included):
```
./main
```
