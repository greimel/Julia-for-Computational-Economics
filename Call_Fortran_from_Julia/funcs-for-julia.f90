module funcs_for_julia
  integer, parameter :: dp = kind(1.d0) ! double precision

contains

!!! The simplest function: take real, return real
  double precision function double_x(x)


    real(dp) :: x

    double_x = 2*x
  end function double_x

!!! Second simplest function: take two reals, change one of them, return real
  double precision function double_x_give_answer(x, answers)

    real(dp) :: x
    real(dp) :: answers(2)

    answers(1) = 42.0
    answers(2) = 42.01

    double_x_give_answer = 2*x
  end function double_x_give_answer

!!! A function that takes a function and returns a value
  double precision function eval_func(f, x)
    real(dp) :: x

    interface
       double precision function f(y)
         integer, parameter :: dp = kind(1.d0) ! double precision
         real(dp) :: y
       end function f
    end interface

    eval_func = f(x)
  end function eval_func

end module funcs_for_julia
