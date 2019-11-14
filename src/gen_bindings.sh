#!/usr/bin/env bash

cat >./fg_interface.h <<EOL
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif
EOL

for n in {6..12}; do
    for table_size in 512 1024 2048 4096 8192; do
        printf "void *fg_init_${n}_${table_size}(double (*)(double), double, double, double, double, uint16_t);\n" >> fg_interface.h
        printf "double fg_eval_${n}_${table_size}(void *, double);\n" >> fg_interface.h
        printf "void fg_delete_${n}_${table_size}(void *);\n\n" >> fg_interface.h
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
EOL

for n in {6..12}; do
    for table_size in 512 1024 2048 4096 8192; do
        printf "typedef FunctionGenerator<${n}, ${table_size}, double> *FGHandle_${n}_${table_size};\n" >> fg_interface.cpp
        printf "void *fg_init_${n}_${table_size}(double (*fin)(double), double a, double b, double tol, double mw, uint16_t error_model) {\n" >> fg_interface.cpp
        printf "    return (void *)new FunctionGenerator<${n}, ${table_size}, double>(fin, a, b, tol, mw, (FGError::ErrorModel) error_model);\n" >> fg_interface.cpp
        printf "}\n" >> fg_interface.cpp

        printf "double fg_eval_${n}_${table_size}(void *f, double x) {\n" >> fg_interface.cpp
        printf "    return (*(FGHandle_${n}_${table_size})f)(x);\n" >> fg_interface.cpp
        printf "}\n" >> fg_interface.cpp
        printf "void fg_delete_${n}_${table_size}(void *f) { delete (FGHandle_${n}_${table_size})f; };\n\n" >> fg_interface.cpp
    done
done

cat <<EOT >> fg_interface.cpp

}
EOT
