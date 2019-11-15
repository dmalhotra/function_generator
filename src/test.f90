
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
  print *, fg_eval(myfun, x)

contains
  function log_wrapper (arg) bind(c) result(y)
    use, intrinsic :: iso_c_binding
    real(kind=c_double), intent(in), value :: arg
    real(kind=c_double) :: y
    y = log(arg)
  end function log_wrapper
  
end program main
