#!/usr/bin/env python3

import numpy as np
from FunctionGenerator import FunctionGenerator as FG
# from function_generator import standard_error_model, relative_error_model
import time
from scipy.special import struve, y0, hankel1
import numba
import os

cpu_count = 1

n = 10000000
approx_range = [1e-10, 1000]
test_range = [1e-10, 999]
tol = 1e-12
order = 12

# functions to test evaluation of
true_funcs = [
    lambda x: y0(x),
    # lambda x: hankel1(0, x),
    lambda x: np.log(x),
    lambda x: 1/x**8,
]
true_func_disp = [
    'y0(x); using standard error model, relative error not guaranteed',
    # 'hankel1(0, x); Using standard error model, relative error not guaranteed',
    'np.log(x); Using standard error model, relative error not guaranteed',
    '1/x**8; Using relative error model, relative error should be good',
]
# error_models = [
#     standard_error_model,
#     standard_error_model,
#     standard_error_model,
#     relative_error_model,
# ]

def random_in(n, a, b):
    x = np.random.rand(n)
    x *= (b-a)
    x += a
    return x

xtest = random_in(n-10000, test_range[0], test_range[1])
xtest = np.concatenate([np.linspace(test_range[0], test_range[1], 10000), xtest])

print('\nTesting function generator on', n, 'points')
print('    minimum test value is: {:0.2e}'.format(xtest.min()))
print('    maximum test value is: {:0.2e}'.format(xtest.max()))
print('')
print('    Standard error model means normalization by max(1, value)')
print('    Relative error model means normalization by value')

for func, disp in zip(true_funcs, true_func_disp):
    print('\n    Function is: ', disp)

    approx_func = FG(func, approx_range[0], approx_range[1], tol, order)
    st = time.time()
    fa = approx_func(xtest)
    approx_func_time1 = time.time()-st

    print('        Approx time (ms):   {:0.1f}'.format(approx_func_time1*1000))
    print('        Points/Sec/Core, Millions:       {:0.1f}'.format(n/approx_func_time1/1000000/cpu_count))

