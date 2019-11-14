#include "fg_interface.h"
#include "function_generator.hpp"

extern "C" {
typedef FunctionGenerator<8, 4096, double> *FGHandle;

void *fg_init_8_2048(double (*fin)(double), double a, double b, double tol,
                     double mw, uint16_t error_model) {
    return (void *)new FunctionGenerator<8, 2048, double>(fin, a, b, tol, mw,
                                                          error_model);
}

void *fg_init_8_4096(double (*fin)(double), double a, double b, double tol,
                     double mw, uint16_t error_model) {
    return (void *)new FunctionGenerator<8, 4096, double>(fin, a, b, tol, mw,
                                                          error_model);
}

double fg_eval(void *f, double x) { return (*(FGHandle)f)(x); }

void fg_delete(void *f) { delete (FGHandle)f; };
}
