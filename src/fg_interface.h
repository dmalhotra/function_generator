#include <inttypes.h>

#ifdef __cplusplus
typedef FunctionGenerator<8, 4096, double> *FGHandle;
extern "C" {
#else
typedef struct FunctionGenerator *FGHandle;
#endif

FGHandle fg_init(double (*)(double), double, double, double, double, uint16_t);
double fg_eval(FGHandle, double);

#ifdef __cplusplus
}
#endif
