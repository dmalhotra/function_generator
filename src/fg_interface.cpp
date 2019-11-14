#include "function_generator.hpp"
#include "fg_interface.h"

extern "C"
{
FGHandle fg_init(double (*fin)(double), double a, double b, double tol,
                     double mw, uint16_t error_model) {
    return new FunctionGenerator<8, 4096, double>(fin, a, b, tol, mw, error_model);
}

double fg_eval(FGHandle f, double x) { return (*f)(x); }
}
