module function_generator
  use, intrinsic :: iso_c_binding
  implicit none

  type, bind(c) :: fg_func
    type(c_ptr) :: obj
    type(c_ptr) :: eval
  end type fg_func

  interface
    function fg_init_8_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_8_4096

    function fg_eval (f, x) result(y)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(fg_func), intent(in) :: f
      real(kind=c_double), intent(in), value :: x
      real(kind=c_double) :: y
    end function fg_eval

  end interface
  
end module function_generator
