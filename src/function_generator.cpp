#include "function_generator.hpp"
#include <complex>

template class FunctionGenerator<8, 4096, double>;
template class FunctionGenerator<8, 4096, std::complex<double>>;

