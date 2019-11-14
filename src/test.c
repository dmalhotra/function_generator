#include "fg_interface.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    fg_func test = fg_init_7_4096(log, 1E-15, 1000, 1E-12, 1E-15, 0);

    int n_el = 100000;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r, 1);

    double *x = (double *)malloc(n_el * sizeof(double));

    for (int i = 0; i < n_el; ++i) {
        x[i] = gsl_rng_uniform(r) * 1000;
    }

    double sum = 0.0;
    for (int i_loop = 0; i_loop < 5000; ++i_loop) {
        for (int i = 0; i < n_el; i++) {
            sum += fg_eval(&test, x[i]);
        }
    }
    printf("%g\n", sum);

    return 0;
}
