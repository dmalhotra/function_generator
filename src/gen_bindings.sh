#!/usr/bin/env bash

N_ARR=(6 7 8 9 10 11 12 13 14)
TABLE_SIZE_ARR=(512 1024 2048 4096 8192)


cat >./fg_interface.h <<EOL
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    void *obj;
    double (*eval)(void*, double);
    void (*free)(void*);
} fg_func;

double fg_eval(fg_func *func, double x);
void fg_free(fg_func *func);

EOL

for n in ${N_ARR[@]}; do
    for table_size in ${TABLE_SIZE_ARR[@]}; do
        printf "fg_func fg_init_${n}_${table_size}(double (*)(double), double, double, double, double, int8_t);\n" >> fg_interface.h
        printf "double fg_eval_${n}_${table_size}(void *, double);\n" >> fg_interface.h
        printf "void fg_free_${n}_${table_size}(void *);\n\n" >> fg_interface.h
    done
done

cat <<EOT >> fg_interface.h
#ifdef __cplusplus
}
#endif
EOT

cat >./fg_interface.cpp <<EOL
#include "fg_interface.h"
#include "function_generator.hpp"

extern "C" {
double fg_eval(fg_func *func, double x) {
    return func->eval(func->obj, x);
};

double fg_eval_(fg_func *func, double x) {
    return func->eval(func->obj, x);
};

void fg_free(fg_func *func) {
    func->free(func->obj);
    func->obj = nullptr;
    func->eval = nullptr;
    func->free = nullptr;
};

void fg_free_(fg_func *func) {
    fg_free(func);
};

EOL

for n in ${N_ARR[@]}; do
    for table_size in ${TABLE_SIZE_ARR[@]}; do
        printf "typedef FunctionGenerator<${n}, ${table_size}, double> *FGHandle_${n}_${table_size};\n" >> fg_interface.cpp
        printf "fg_func fg_init_${n}_${table_size}(double (*fin)(double), double a, double b, double tol, double mw, int8_t error_model) {\n" >> fg_interface.cpp
        printf "    fg_func res;\n" >> fg_interface.cpp
        printf "    res.obj = (void *)new FunctionGenerator<${n}, ${table_size}, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);\n" >> fg_interface.cpp
        printf "    res.eval = &fg_eval_${n}_${table_size};\n" >> fg_interface.cpp
        printf "    res.free = &fg_free_${n}_${table_size};\n" >> fg_interface.cpp
        printf "    return res;\n" >> fg_interface.cpp
        printf "}\n" >> fg_interface.cpp

        printf "double fg_eval_${n}_${table_size}(void *f, double x) {\n" >> fg_interface.cpp
        printf "    return (*(FGHandle_${n}_${table_size})f)(x);\n" >> fg_interface.cpp
        printf "}\n" >> fg_interface.cpp
        printf "void fg_free_${n}_${table_size}(void *f) { delete (FGHandle_${n}_${table_size})f; };\n\n" >> fg_interface.cpp
    done
done

cat <<EOT >> fg_interface.cpp

}
EOT

cat >./fg_interface.f90 <<EOL
module function_generator
  use, intrinsic :: iso_c_binding
  implicit none

  type, bind(c) :: fg_func
    type(c_ptr) :: obj
    type(c_ptr) :: eval
    type(c_ptr) :: free
  end type fg_func

  interface
EOL

for n in ${N_ARR[@]}; do
    for table_size in ${TABLE_SIZE_ARR[@]}; do
        printf "    function fg_init_${n}_${table_size} (fin, a, b, tol, mw, error_model) bind(c) result(func)\n" >> fg_interface.f90
        printf "      use, intrinsic :: iso_c_binding\n" >> fg_interface.f90
        printf "      import fg_func\n" >> fg_interface.f90
        printf "      type(c_funptr), intent(in), value :: fin\n" >> fg_interface.f90
        printf "      real(kind=c_double), intent(in), value :: a, b, tol, mw\n" >> fg_interface.f90
        printf "      integer(kind=c_int8_t), intent(in), value :: error_model\n" >> fg_interface.f90
        printf "      type(fg_func) :: func\n" >> fg_interface.f90
        printf "    end function fg_init_${n}_${table_size}\n" >> fg_interface.f90
        printf "\n" >> fg_interface.f90
    done
done

cat <<EOT >> fg_interface.f90
    function fg_eval (f, x) result(y)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(fg_func), intent(in) :: f
      real(kind=c_double), intent(in), value :: x
      real(kind=c_double) :: y
    end function fg_eval
  end interface

end module function_generator
EOT
