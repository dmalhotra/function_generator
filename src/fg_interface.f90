module function_generator
  use, intrinsic :: iso_c_binding
  implicit none

  type, bind(c) :: fg_func
    type(c_ptr) :: obj
    type(c_ptr) :: eval
  end type fg_func

  interface
    function fg_init_6_512 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_6_512

    function fg_init_6_1024 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_6_1024

    function fg_init_6_2048 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_6_2048

    function fg_init_6_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_6_4096

    function fg_init_6_8192 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_6_8192

    function fg_init_7_512 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_7_512

    function fg_init_7_1024 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_7_1024

    function fg_init_7_2048 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_7_2048

    function fg_init_7_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_7_4096

    function fg_init_7_8192 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_7_8192

    function fg_init_8_512 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_8_512

    function fg_init_8_1024 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_8_1024

    function fg_init_8_2048 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_8_2048

    function fg_init_8_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_8_4096

    function fg_init_8_8192 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_8_8192

    function fg_init_9_512 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_9_512

    function fg_init_9_1024 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_9_1024

    function fg_init_9_2048 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_9_2048

    function fg_init_9_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_9_4096

    function fg_init_9_8192 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_9_8192

    function fg_init_10_512 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_10_512

    function fg_init_10_1024 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_10_1024

    function fg_init_10_2048 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_10_2048

    function fg_init_10_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_10_4096

    function fg_init_10_8192 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_10_8192

    function fg_init_11_512 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_11_512

    function fg_init_11_1024 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_11_1024

    function fg_init_11_2048 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_11_2048

    function fg_init_11_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_11_4096

    function fg_init_11_8192 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_11_8192

    function fg_init_12_512 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_12_512

    function fg_init_12_1024 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_12_1024

    function fg_init_12_2048 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_12_2048

    function fg_init_12_4096 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_12_4096

    function fg_init_12_8192 (fin, a, b, tol, mw, error_model) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(c_funptr), intent(in), value :: fin
      real(kind=c_double), intent(in), value :: a, b, tol, mw
      integer(kind=c_int8_t), intent(in), value :: error_model
      type(fg_func) :: func
    end function fg_init_12_8192

    function fg_eval (f, x) result(y)
      use, intrinsic :: iso_c_binding
      import fg_func
      type(fg_func), intent(in) :: f
      real(kind=c_double), intent(in), value :: x
      real(kind=c_double) :: y
    end function fg_eval

  end interface
  
end module function_generator
