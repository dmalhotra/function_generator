#include "fg_interface.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    FGHandle f = fg_init(log, 1E-15, 1000, 1E-12, 1E-15, 0);
    printf("%g\n", fg_eval(f, 0.5));
}
