#include "fg_interface.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    void *f1 = fg_init_8_4096(log, 1E-15, 1000, 1E-12, 1E-15, 0);
    void *f2 = fg_init_8_2048(log, 1E-15, 1000, 1E-12, 1E-15, 0);

    printf("%g\n", fg_eval(f1, 0.5));
    printf("%g\n", fg_eval(f2, 0.5));

    fg_delete(f1);
    fg_delete(f2);
    return 0;
}
