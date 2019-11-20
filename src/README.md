# Experimental C++ variant of FunctionGenerator with C/Fortran bindings

## General warning
First, this project is a work in progress, so the API and functionality is likely to change. If
you have any suggestions for functionality, performance, API, or general useability changes,
please don't hesitate to include a pull request, create an issue, or contact me directly. Also,
if you find this project useful at all, please let me know.

## About the project
I originally developed this project to speed up @dbstein's excellent `python` FunctionGenerator
package, which this project has been forked from. This fork resulted in some algorithmic
improvements to both of our packages, but ultimately the statically generated `C++` code could
not typically compete with `numba`. This is most likely due to the `llvm` `jit` compiler in
numba having more information at compile time than can be provided in a statically compiled
language such as `C++`. This is almost entirely alleviated by the, at first glance, bizarre
choice to template some of the parameters of the class (expansion order and size of lookup
table). This choice reduces the efficacy of python bindings by having to explicitly bind every
variant of the class to a unique name (an issue that is shared by `fortran` and `C`). Given
that the performance of the `C++` library is, at best, on par with the python variant, the
python bindings were removed and instead the library is developed separately but with added
`fortran` and `C` bindings.

## How it works
The algorithm is straight forward. During initialization, the `FunctionGenerator` attempts to
fit the target function to a Chebyshev polynomial to a given input order. If it is unable to
fit the target function to a given tolerance, the region is split evenly into two sub-divisions
and the process repeats for each sub-division. This process is repeated recursively until every
sub-division meets the given error tolerance.

When the resultant approximation needs to be evaluated at say, point `x`, the algorithm needs
to quickly find the correct sub-division to use for the approximation. Given that there can be
hundreds of sub-divisions which are extremely heterogenous in size, a simple bisection routine
can be quite slow, taking sometimes much more time than the actual expansion calculation. The
algorithm is therefore aided by a lookup table to help quickly find a smaller set of
sub-divisions that `x` could lie within. For most functions across most values, this process
will typically be `O(1)`, though still in worst case `N log(N)`, where `N` is the number of
sub-divisions.

## Using the package
### Parameters
The version as of writing this documentation has eight inputs.
1. `n`, the order of the Chebyshev expansion. Low orders evaluate faster, but will typically
   require more subdivisions and is more likely miss features of the underlying function due to
   how the error is currently measured. This typically only happens if there is a spike in the
   function far from any Chebyshev nodes. This will likely change in the future to be more
   robust.  Values of 6-12 have proven most effective for the given test functions, with lower
   values more often outperforming high values. This is a template parameter in the `C++`
   class, and varies linearly from 6 to 14 in the `C/fortran` bindings
2. `table_size`, the number of elements in the lookup table to speed up the lookup
   process. This is a template parameter in the `C++` class, and is compiled for values in
   powers of 2 from 512 to 8192 in the `C/fortran` bindings.
3. `fin`, an input function. A one dimension function that returns a `double` and takes a
   `double` as its sole input argument. It is templated in the `C++` version to return
   non-floating point values, but this is currently broken/unsupported. If it is needed, it
   should be a relatively quick fix to include, but just provides extra pain in maintaining the
   bindings. See below codes for examples of how to provide input functions.
4. `low`, the lower bound of your input domain.
5. `high`, the upper bound of your input domain.
6. `tol`, the desired accuracy of your function. See main `README` of parent project (directory
   up from the current directory) for more details.
7. `minimum_width`, the minimum width for an allowed sub-division. If this is exceeded during
   initialization, the program should abort and throw an error message.
8. `error_model`, error model to use to judge if within `tol` parameter. Currently only
   standard (0), and relative (0) are supported. See `README` of parent project for more on
   this.

### Dependencies
The only dependency is `Eigen3`. The `C` example also requires `gsl`. If eigen is not in the
default include path for your compiler, you must supply it. For the examples just use, assuming
`EIGEN_BASE` is defined as the path to your eigen installation, `CFLAGS=${EIGEN_BASE}/include
make all` should work. If you have installed these dependencies via `conda`, then the make file
should automatically use these versions.

### C++
In C++, no compilation is necessary, though pre-compilation into a shared library is also
possible if desired. At least `c++11` is required.
```
// test.cpp
#include <iostream>
#include "function_generator.hpp"

int main(int argc, char *argv[]) {
    constexpr int n = 8;
    constexpr int table_size = 4096;
    auto fin = static_cast<double (*)(double)>(std::log); // Function to interpolate
    double low = 1E-15; // Lower bound of function's domain to interpolate
    double high = 1000; // Upper bound of function's domain to interpolate
    double tol = 1E-10; // Desired accuracy
    double minimum_width = 1E-15; // Minimum width of interpolation sub-region
    auto error_model = FGError::ErrorModel::standard;

    FunctionGenerator<n, table_size, double> f(fin, low, high, tol,
                                               minimum_width, error_model);

    std::cout << std::fabs(f(1.5) - fin(1.5)) << std::endl;
    return 0;
}
```

```
c++ test.cpp -I${EIGEN_BASE}/include -std=c++11
./a.out
>>> 6.05072e-15
```

### C
Using the `C` bindings requires compilation of an intermediate 'interface' object defined by
`fg_interface.cpp` and `fg_interface.h`. You can see the `Makefile` for an example, but a more
basic one is outlined here as well.

There are a few things to note first.
1. There are a *lot* of functions. This is because, since the library is templated, there must
   be a separate class compiled for each used template. The resulting compiled objects will
   possibly be very different from one another, and so can not be used interchangeably. The
   function naming convention is `fg_${operation}_${expansion_order}_${table_size}`. When the
   function needs the underlying object (everything but `init`), the first argument should
   always be a pointer to the struct returned by the `init` routine.
2. There is a wrapper struct that helps to mimic the object oriented `C++` class I have named
   `fg_func`. It contains a pointer to the `C++` object `fg_func.obj`, the relevant evaluation
   function `fg_func.eval`, and the relevant deletion routine `fg_func.free`. This struct is
   subject to have more data in the future.
3. I have provided the `fg_eval` convenience wrapper. This has a mild function pointer lookup
   overhead compared to direct calls like `fg_eval_8_4096(test.obj, x)`, which can be used
   instead for the best performance. In my tests, the results were typically a few percent
   faster using the direct function call. There should be no issue using the "fully qualified"
   bindings, but if you call a different function than the one intended by the init call, your
   results will very likely be wrong, or the program might crash, if you're lucky. Use with caution!
4. I have also provided the `fg_free` convenience wrapper. There really is no reason to use the
   `fg_free_${n}_${ts}` bindings directly, but they are there.

```
// test.c
#include "fg_interface.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
    // Function naming convention (for doubles) is fg_init_${expansion_order}_${table_size} and then
    // the standard as defined in function_generator.hpp
    fg_func test = fg_init_8_4096(log, 1E-15, 1000, 1E-12, 1E-15, 0);

    printf("%g\n", fabs(fg_eval(&test, 1.5) - log(1.5)));

    fg_free(&test);

    return 0;
}
```

```
c++ -std=c++11 -c fg_interface.cpp
cc -c test.c
cc test.o fg_interface.o -lstdc++ -lm
./a.out >>> 0
```

### Fortran
I don't know `fortran` beyond the very basics, but I have provided some `fortran` bindings to
the `C` bindings the best of my ability. Since the two sets of bindings are basically the same,
you can see the `C` section for calling details.

```
! test.f90
program main
  use function_generator
  implicit none

  type(c_funptr) :: cproc
  real(kind=c_double) :: a, b, tol, mw, x
  integer(kind=c_int8_t) :: error_model
  type(fg_func) :: myfun

  cproc = c_funloc(log_wrapper)
  a = 1E-15
  b = 1000
  tol = 1E-10
  mw = 1E-15
  error_model = 0

  myfun = fg_init_8_4096(cproc, a, b, tol, mw, error_model)
  x = 0.5
  print *, fg_eval(myfun, x) - log(x)

  call fg_free(myfun)

contains
  function log_wrapper (arg) bind(c) result(y)
    use, intrinsic :: iso_c_binding
    real(kind=c_double), intent(in), value :: arg
    real(kind=c_double) :: y
    y = log(arg)
  end function log_wrapper

end program main
```

*NOTE* if you compile `fg_interface.f90` and `fg_interface.cpp`, with the standard `-c` option
to generate objects with no specified output file, the `.o` files will collide. This is why I
explicitly output the compiled `fg_interface.cpp` file to `function_generator.o`. This will
likely change in the future because it's just asking for trouble.

```
g++ -std=c++11 -c fg_interface.cpp -o function_generator.o
f95 -c fg_interface.f90
f95 -c test.f90
f95 function_generator.o fg_interface.o test.o -lstdc++
./a.out
>>>  -3.3306690738754696E-016
```
