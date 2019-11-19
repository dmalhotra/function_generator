#include "fg_interface.h"
#include "function_generator.hpp"

extern "C" {
// Fortran binding
double fg_eval_(fg_func *func, double x) {
    return func->eval(func->obj, x);
};

typedef FunctionGenerator<6, 512, double> *FGHandle_6_512;
fg_func fg_init_6_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<6, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_6_512;
    return res;
}
double fg_eval_6_512(void *f, double x) {
    return (*(FGHandle_6_512)f)(x);
}
void fg_delete_6_512(void *f) { delete (FGHandle_6_512)f; };

typedef FunctionGenerator<6, 1024, double> *FGHandle_6_1024;
fg_func fg_init_6_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<6, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_6_1024;
    return res;
}
double fg_eval_6_1024(void *f, double x) {
    return (*(FGHandle_6_1024)f)(x);
}
void fg_delete_6_1024(void *f) { delete (FGHandle_6_1024)f; };

typedef FunctionGenerator<6, 2048, double> *FGHandle_6_2048;
fg_func fg_init_6_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<6, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_6_2048;
    return res;
}
double fg_eval_6_2048(void *f, double x) {
    return (*(FGHandle_6_2048)f)(x);
}
void fg_delete_6_2048(void *f) { delete (FGHandle_6_2048)f; };

typedef FunctionGenerator<6, 4096, double> *FGHandle_6_4096;
fg_func fg_init_6_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<6, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_6_4096;
    return res;
}
double fg_eval_6_4096(void *f, double x) {
    return (*(FGHandle_6_4096)f)(x);
}
void fg_delete_6_4096(void *f) { delete (FGHandle_6_4096)f; };

typedef FunctionGenerator<6, 8192, double> *FGHandle_6_8192;
fg_func fg_init_6_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<6, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_6_8192;
    return res;
}
double fg_eval_6_8192(void *f, double x) {
    return (*(FGHandle_6_8192)f)(x);
}
void fg_delete_6_8192(void *f) { delete (FGHandle_6_8192)f; };

typedef FunctionGenerator<7, 512, double> *FGHandle_7_512;
fg_func fg_init_7_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<7, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_7_512;
    return res;
}
double fg_eval_7_512(void *f, double x) {
    return (*(FGHandle_7_512)f)(x);
}
void fg_delete_7_512(void *f) { delete (FGHandle_7_512)f; };

typedef FunctionGenerator<7, 1024, double> *FGHandle_7_1024;
fg_func fg_init_7_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<7, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_7_1024;
    return res;
}
double fg_eval_7_1024(void *f, double x) {
    return (*(FGHandle_7_1024)f)(x);
}
void fg_delete_7_1024(void *f) { delete (FGHandle_7_1024)f; };

typedef FunctionGenerator<7, 2048, double> *FGHandle_7_2048;
fg_func fg_init_7_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<7, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_7_2048;
    return res;
}
double fg_eval_7_2048(void *f, double x) {
    return (*(FGHandle_7_2048)f)(x);
}
void fg_delete_7_2048(void *f) { delete (FGHandle_7_2048)f; };

typedef FunctionGenerator<7, 4096, double> *FGHandle_7_4096;
fg_func fg_init_7_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<7, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_7_4096;
    return res;
}
double fg_eval_7_4096(void *f, double x) {
    return (*(FGHandle_7_4096)f)(x);
}
void fg_delete_7_4096(void *f) { delete (FGHandle_7_4096)f; };

typedef FunctionGenerator<7, 8192, double> *FGHandle_7_8192;
fg_func fg_init_7_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<7, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_7_8192;
    return res;
}
double fg_eval_7_8192(void *f, double x) {
    return (*(FGHandle_7_8192)f)(x);
}
void fg_delete_7_8192(void *f) { delete (FGHandle_7_8192)f; };

typedef FunctionGenerator<8, 512, double> *FGHandle_8_512;
fg_func fg_init_8_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<8, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_8_512;
    return res;
}
double fg_eval_8_512(void *f, double x) {
    return (*(FGHandle_8_512)f)(x);
}
void fg_delete_8_512(void *f) { delete (FGHandle_8_512)f; };

typedef FunctionGenerator<8, 1024, double> *FGHandle_8_1024;
fg_func fg_init_8_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<8, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_8_1024;
    return res;
}
double fg_eval_8_1024(void *f, double x) {
    return (*(FGHandle_8_1024)f)(x);
}
void fg_delete_8_1024(void *f) { delete (FGHandle_8_1024)f; };

typedef FunctionGenerator<8, 2048, double> *FGHandle_8_2048;
fg_func fg_init_8_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<8, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_8_2048;
    return res;
}
double fg_eval_8_2048(void *f, double x) {
    return (*(FGHandle_8_2048)f)(x);
}
void fg_delete_8_2048(void *f) { delete (FGHandle_8_2048)f; };

typedef FunctionGenerator<8, 4096, double> *FGHandle_8_4096;
fg_func fg_init_8_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<8, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_8_4096;
    return res;
}
double fg_eval_8_4096(void *f, double x) {
    return (*(FGHandle_8_4096)f)(x);
}
void fg_delete_8_4096(void *f) { delete (FGHandle_8_4096)f; };

typedef FunctionGenerator<8, 8192, double> *FGHandle_8_8192;
fg_func fg_init_8_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<8, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_8_8192;
    return res;
}
double fg_eval_8_8192(void *f, double x) {
    return (*(FGHandle_8_8192)f)(x);
}
void fg_delete_8_8192(void *f) { delete (FGHandle_8_8192)f; };

typedef FunctionGenerator<9, 512, double> *FGHandle_9_512;
fg_func fg_init_9_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<9, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_9_512;
    return res;
}
double fg_eval_9_512(void *f, double x) {
    return (*(FGHandle_9_512)f)(x);
}
void fg_delete_9_512(void *f) { delete (FGHandle_9_512)f; };

typedef FunctionGenerator<9, 1024, double> *FGHandle_9_1024;
fg_func fg_init_9_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<9, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_9_1024;
    return res;
}
double fg_eval_9_1024(void *f, double x) {
    return (*(FGHandle_9_1024)f)(x);
}
void fg_delete_9_1024(void *f) { delete (FGHandle_9_1024)f; };

typedef FunctionGenerator<9, 2048, double> *FGHandle_9_2048;
fg_func fg_init_9_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<9, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_9_2048;
    return res;
}
double fg_eval_9_2048(void *f, double x) {
    return (*(FGHandle_9_2048)f)(x);
}
void fg_delete_9_2048(void *f) { delete (FGHandle_9_2048)f; };

typedef FunctionGenerator<9, 4096, double> *FGHandle_9_4096;
fg_func fg_init_9_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<9, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_9_4096;
    return res;
}
double fg_eval_9_4096(void *f, double x) {
    return (*(FGHandle_9_4096)f)(x);
}
void fg_delete_9_4096(void *f) { delete (FGHandle_9_4096)f; };

typedef FunctionGenerator<9, 8192, double> *FGHandle_9_8192;
fg_func fg_init_9_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<9, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_9_8192;
    return res;
}
double fg_eval_9_8192(void *f, double x) {
    return (*(FGHandle_9_8192)f)(x);
}
void fg_delete_9_8192(void *f) { delete (FGHandle_9_8192)f; };

typedef FunctionGenerator<10, 512, double> *FGHandle_10_512;
fg_func fg_init_10_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<10, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_10_512;
    return res;
}
double fg_eval_10_512(void *f, double x) {
    return (*(FGHandle_10_512)f)(x);
}
void fg_delete_10_512(void *f) { delete (FGHandle_10_512)f; };

typedef FunctionGenerator<10, 1024, double> *FGHandle_10_1024;
fg_func fg_init_10_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<10, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_10_1024;
    return res;
}
double fg_eval_10_1024(void *f, double x) {
    return (*(FGHandle_10_1024)f)(x);
}
void fg_delete_10_1024(void *f) { delete (FGHandle_10_1024)f; };

typedef FunctionGenerator<10, 2048, double> *FGHandle_10_2048;
fg_func fg_init_10_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<10, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_10_2048;
    return res;
}
double fg_eval_10_2048(void *f, double x) {
    return (*(FGHandle_10_2048)f)(x);
}
void fg_delete_10_2048(void *f) { delete (FGHandle_10_2048)f; };

typedef FunctionGenerator<10, 4096, double> *FGHandle_10_4096;
fg_func fg_init_10_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<10, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_10_4096;
    return res;
}
double fg_eval_10_4096(void *f, double x) {
    return (*(FGHandle_10_4096)f)(x);
}
void fg_delete_10_4096(void *f) { delete (FGHandle_10_4096)f; };

typedef FunctionGenerator<10, 8192, double> *FGHandle_10_8192;
fg_func fg_init_10_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<10, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_10_8192;
    return res;
}
double fg_eval_10_8192(void *f, double x) {
    return (*(FGHandle_10_8192)f)(x);
}
void fg_delete_10_8192(void *f) { delete (FGHandle_10_8192)f; };

typedef FunctionGenerator<11, 512, double> *FGHandle_11_512;
fg_func fg_init_11_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<11, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_11_512;
    return res;
}
double fg_eval_11_512(void *f, double x) {
    return (*(FGHandle_11_512)f)(x);
}
void fg_delete_11_512(void *f) { delete (FGHandle_11_512)f; };

typedef FunctionGenerator<11, 1024, double> *FGHandle_11_1024;
fg_func fg_init_11_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<11, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_11_1024;
    return res;
}
double fg_eval_11_1024(void *f, double x) {
    return (*(FGHandle_11_1024)f)(x);
}
void fg_delete_11_1024(void *f) { delete (FGHandle_11_1024)f; };

typedef FunctionGenerator<11, 2048, double> *FGHandle_11_2048;
fg_func fg_init_11_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<11, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_11_2048;
    return res;
}
double fg_eval_11_2048(void *f, double x) {
    return (*(FGHandle_11_2048)f)(x);
}
void fg_delete_11_2048(void *f) { delete (FGHandle_11_2048)f; };

typedef FunctionGenerator<11, 4096, double> *FGHandle_11_4096;
fg_func fg_init_11_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<11, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_11_4096;
    return res;
}
double fg_eval_11_4096(void *f, double x) {
    return (*(FGHandle_11_4096)f)(x);
}
void fg_delete_11_4096(void *f) { delete (FGHandle_11_4096)f; };

typedef FunctionGenerator<11, 8192, double> *FGHandle_11_8192;
fg_func fg_init_11_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<11, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_11_8192;
    return res;
}
double fg_eval_11_8192(void *f, double x) {
    return (*(FGHandle_11_8192)f)(x);
}
void fg_delete_11_8192(void *f) { delete (FGHandle_11_8192)f; };

typedef FunctionGenerator<12, 512, double> *FGHandle_12_512;
fg_func fg_init_12_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<12, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_12_512;
    return res;
}
double fg_eval_12_512(void *f, double x) {
    return (*(FGHandle_12_512)f)(x);
}
void fg_delete_12_512(void *f) { delete (FGHandle_12_512)f; };

typedef FunctionGenerator<12, 1024, double> *FGHandle_12_1024;
fg_func fg_init_12_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<12, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_12_1024;
    return res;
}
double fg_eval_12_1024(void *f, double x) {
    return (*(FGHandle_12_1024)f)(x);
}
void fg_delete_12_1024(void *f) { delete (FGHandle_12_1024)f; };

typedef FunctionGenerator<12, 2048, double> *FGHandle_12_2048;
fg_func fg_init_12_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<12, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_12_2048;
    return res;
}
double fg_eval_12_2048(void *f, double x) {
    return (*(FGHandle_12_2048)f)(x);
}
void fg_delete_12_2048(void *f) { delete (FGHandle_12_2048)f; };

typedef FunctionGenerator<12, 4096, double> *FGHandle_12_4096;
fg_func fg_init_12_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<12, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_12_4096;
    return res;
}
double fg_eval_12_4096(void *f, double x) {
    return (*(FGHandle_12_4096)f)(x);
}
void fg_delete_12_4096(void *f) { delete (FGHandle_12_4096)f; };

typedef FunctionGenerator<12, 8192, double> *FGHandle_12_8192;
fg_func fg_init_12_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<12, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_12_8192;
    return res;
}
double fg_eval_12_8192(void *f, double x) {
    return (*(FGHandle_12_8192)f)(x);
}
void fg_delete_12_8192(void *f) { delete (FGHandle_12_8192)f; };

typedef FunctionGenerator<13, 512, double> *FGHandle_13_512;
fg_func fg_init_13_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<13, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_13_512;
    return res;
}
double fg_eval_13_512(void *f, double x) {
    return (*(FGHandle_13_512)f)(x);
}
void fg_delete_13_512(void *f) { delete (FGHandle_13_512)f; };

typedef FunctionGenerator<13, 1024, double> *FGHandle_13_1024;
fg_func fg_init_13_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<13, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_13_1024;
    return res;
}
double fg_eval_13_1024(void *f, double x) {
    return (*(FGHandle_13_1024)f)(x);
}
void fg_delete_13_1024(void *f) { delete (FGHandle_13_1024)f; };

typedef FunctionGenerator<13, 2048, double> *FGHandle_13_2048;
fg_func fg_init_13_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<13, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_13_2048;
    return res;
}
double fg_eval_13_2048(void *f, double x) {
    return (*(FGHandle_13_2048)f)(x);
}
void fg_delete_13_2048(void *f) { delete (FGHandle_13_2048)f; };

typedef FunctionGenerator<13, 4096, double> *FGHandle_13_4096;
fg_func fg_init_13_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<13, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_13_4096;
    return res;
}
double fg_eval_13_4096(void *f, double x) {
    return (*(FGHandle_13_4096)f)(x);
}
void fg_delete_13_4096(void *f) { delete (FGHandle_13_4096)f; };

typedef FunctionGenerator<13, 8192, double> *FGHandle_13_8192;
fg_func fg_init_13_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<13, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_13_8192;
    return res;
}
double fg_eval_13_8192(void *f, double x) {
    return (*(FGHandle_13_8192)f)(x);
}
void fg_delete_13_8192(void *f) { delete (FGHandle_13_8192)f; };

typedef FunctionGenerator<14, 512, double> *FGHandle_14_512;
fg_func fg_init_14_512(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<14, 512, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_14_512;
    return res;
}
double fg_eval_14_512(void *f, double x) {
    return (*(FGHandle_14_512)f)(x);
}
void fg_delete_14_512(void *f) { delete (FGHandle_14_512)f; };

typedef FunctionGenerator<14, 1024, double> *FGHandle_14_1024;
fg_func fg_init_14_1024(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<14, 1024, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_14_1024;
    return res;
}
double fg_eval_14_1024(void *f, double x) {
    return (*(FGHandle_14_1024)f)(x);
}
void fg_delete_14_1024(void *f) { delete (FGHandle_14_1024)f; };

typedef FunctionGenerator<14, 2048, double> *FGHandle_14_2048;
fg_func fg_init_14_2048(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<14, 2048, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_14_2048;
    return res;
}
double fg_eval_14_2048(void *f, double x) {
    return (*(FGHandle_14_2048)f)(x);
}
void fg_delete_14_2048(void *f) { delete (FGHandle_14_2048)f; };

typedef FunctionGenerator<14, 4096, double> *FGHandle_14_4096;
fg_func fg_init_14_4096(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<14, 4096, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_14_4096;
    return res;
}
double fg_eval_14_4096(void *f, double x) {
    return (*(FGHandle_14_4096)f)(x);
}
void fg_delete_14_4096(void *f) { delete (FGHandle_14_4096)f; };

typedef FunctionGenerator<14, 8192, double> *FGHandle_14_8192;
fg_func fg_init_14_8192(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {
    fg_func res;
    res.obj = (void *)new FunctionGenerator<14, 8192, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);
    res.eval = &fg_eval_14_8192;
    return res;
}
double fg_eval_14_8192(void *f, double x) {
    return (*(FGHandle_14_8192)f)(x);
}
void fg_delete_14_8192(void *f) { delete (FGHandle_14_8192)f; };


}
