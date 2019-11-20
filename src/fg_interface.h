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

fg_func fg_init_6_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_6_512(void *, double);
void fg_free_6_512(void *);

fg_func fg_init_6_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_6_1024(void *, double);
void fg_free_6_1024(void *);

fg_func fg_init_6_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_6_2048(void *, double);
void fg_free_6_2048(void *);

fg_func fg_init_6_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_6_4096(void *, double);
void fg_free_6_4096(void *);

fg_func fg_init_6_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_6_8192(void *, double);
void fg_free_6_8192(void *);

fg_func fg_init_7_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_7_512(void *, double);
void fg_free_7_512(void *);

fg_func fg_init_7_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_7_1024(void *, double);
void fg_free_7_1024(void *);

fg_func fg_init_7_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_7_2048(void *, double);
void fg_free_7_2048(void *);

fg_func fg_init_7_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_7_4096(void *, double);
void fg_free_7_4096(void *);

fg_func fg_init_7_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_7_8192(void *, double);
void fg_free_7_8192(void *);

fg_func fg_init_8_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_8_512(void *, double);
void fg_free_8_512(void *);

fg_func fg_init_8_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_8_1024(void *, double);
void fg_free_8_1024(void *);

fg_func fg_init_8_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_8_2048(void *, double);
void fg_free_8_2048(void *);

fg_func fg_init_8_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_8_4096(void *, double);
void fg_free_8_4096(void *);

fg_func fg_init_8_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_8_8192(void *, double);
void fg_free_8_8192(void *);

fg_func fg_init_9_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_9_512(void *, double);
void fg_free_9_512(void *);

fg_func fg_init_9_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_9_1024(void *, double);
void fg_free_9_1024(void *);

fg_func fg_init_9_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_9_2048(void *, double);
void fg_free_9_2048(void *);

fg_func fg_init_9_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_9_4096(void *, double);
void fg_free_9_4096(void *);

fg_func fg_init_9_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_9_8192(void *, double);
void fg_free_9_8192(void *);

fg_func fg_init_10_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_10_512(void *, double);
void fg_free_10_512(void *);

fg_func fg_init_10_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_10_1024(void *, double);
void fg_free_10_1024(void *);

fg_func fg_init_10_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_10_2048(void *, double);
void fg_free_10_2048(void *);

fg_func fg_init_10_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_10_4096(void *, double);
void fg_free_10_4096(void *);

fg_func fg_init_10_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_10_8192(void *, double);
void fg_free_10_8192(void *);

fg_func fg_init_11_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_11_512(void *, double);
void fg_free_11_512(void *);

fg_func fg_init_11_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_11_1024(void *, double);
void fg_free_11_1024(void *);

fg_func fg_init_11_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_11_2048(void *, double);
void fg_free_11_2048(void *);

fg_func fg_init_11_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_11_4096(void *, double);
void fg_free_11_4096(void *);

fg_func fg_init_11_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_11_8192(void *, double);
void fg_free_11_8192(void *);

fg_func fg_init_12_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_12_512(void *, double);
void fg_free_12_512(void *);

fg_func fg_init_12_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_12_1024(void *, double);
void fg_free_12_1024(void *);

fg_func fg_init_12_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_12_2048(void *, double);
void fg_free_12_2048(void *);

fg_func fg_init_12_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_12_4096(void *, double);
void fg_free_12_4096(void *);

fg_func fg_init_12_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_12_8192(void *, double);
void fg_free_12_8192(void *);

fg_func fg_init_13_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_13_512(void *, double);
void fg_free_13_512(void *);

fg_func fg_init_13_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_13_1024(void *, double);
void fg_free_13_1024(void *);

fg_func fg_init_13_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_13_2048(void *, double);
void fg_free_13_2048(void *);

fg_func fg_init_13_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_13_4096(void *, double);
void fg_free_13_4096(void *);

fg_func fg_init_13_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_13_8192(void *, double);
void fg_free_13_8192(void *);

fg_func fg_init_14_512(double (*)(double), double, double, double, double, int8_t);
double fg_eval_14_512(void *, double);
void fg_free_14_512(void *);

fg_func fg_init_14_1024(double (*)(double), double, double, double, double, int8_t);
double fg_eval_14_1024(void *, double);
void fg_free_14_1024(void *);

fg_func fg_init_14_2048(double (*)(double), double, double, double, double, int8_t);
double fg_eval_14_2048(void *, double);
void fg_free_14_2048(void *);

fg_func fg_init_14_4096(double (*)(double), double, double, double, double, int8_t);
double fg_eval_14_4096(void *, double);
void fg_free_14_4096(void *);

fg_func fg_init_14_8192(double (*)(double), double, double, double, double, int8_t);
double fg_eval_14_8192(void *, double);
void fg_free_14_8192(void *);

#ifdef __cplusplus
}
#endif
