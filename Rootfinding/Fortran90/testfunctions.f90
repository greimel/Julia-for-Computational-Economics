! module counter
!   use constants
!   integer(i4b) :: neval
! end module counter

module functions
  use constants
!  use counter
  
contains
  double precision function f1(x)
    use constants
    real(ndp) :: a, b, c, d, x
    a = 0.5
    b = 0.5
    c = -0.5
    d = -0.2

    ! neval = neval + 1
    f1 = a * x**c + b*x**d

    return
  end function f1

  double precision function f1p(x)
    use constants
    real(ndp) :: a, b, c, d, x
    a = 0.5
    b = 0.5
    c = -0.5
    d = -0.2

    ! neval = neval + 1
    f1p = a * c * x**(c-1.0) + b*d*x**(d-1)

    return
  end function f1p

  double precision function f2(x)
    use constants
    real(ndp) :: a, b, x
    a = 1.5
    b = 1.0

   ! neval = neval + 1
    f2 = dlog(a + x) - b
    
    return
  end function f2

  double precision function f2p(x)
    use constants
    real(ndp) :: a, b, x
    a = 1.5
    b = 1.0
   ! neval = neval + 1
    f2p = 1.0 / (a + x)

    return
  end function f2p

  double precision function f3(x)
    use constants
    real(ndp) :: x, a0, y0, y1, beta, r, sigma
    a0 = 1
    y0 = 1
    y1 = 1
    beta = 0.99
    r = 1.03
    sigma = 2.0

    !neval = neval + 1

    f3 = - (a0 + y0 - x) ** (-sigma) + beta * (1+r) * ((1+r)*x + y1) ** (-sigma)
    return
  end function f3

  
  double precision function f4(x)
    use constants
  !  use lambdamod
    real(ndp) :: x, a0, y0, ybar, yH, yL, beta, r, sigma, lam
    a0 = 1
    y0 = 1
    ybar = 1
    beta = 0.99
    r = 1.03
    sigma = 2.0
    lam = 0.25

    yL = ybar * (1 - lam)
    yH = ybar * (1 + lam)
    
    !neval = neval + 1

    f4 = - (a0 + y0 - x) ** (-sigma) + beta * (1+r) / 2.0 * (((1+r)*x + yL) ** (-sigma) + ((1+r)*x + yH) ** (-sigma))
    return
  end function f4
  
end module functions
