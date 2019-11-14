#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

void *fg_init_8_2048(double (*)(double), double, double, double, double,
                     uint16_t);
void *fg_init_8_4096(double (*)(double), double, double, double, double,
                     uint16_t);

double fg_eval(void *f, double x);
void fg_delete(void *f);
#ifdef __cplusplus
}
#endif
