! Solve a nonlinear equation using bisection.
! written by Fabian Greimel
! Adapted from Tony Smith

module rootfinding

  use constants

contains

  double precision function bisect(f, xl, xh, mxiter, toler, verbose)

    use constants
    real(ndp) :: xl, xh, xlow, xhigh, xcur, tmp
    real(ndp), intent(in) :: toler
    integer(i4b), intent(in) :: mxiter
    logical :: verbose
    
    interface
       double precision function f(x)
         use constants
         real(ndp) :: x
       end function f
    end interface

    xlow = xl
    xhigh = xh
    
    !write(6,"('start ... bracket: ',3f15.8)") xlow,xhigh,f((xlow+xhigh)/2.0)

    if (f(xcur) > 0) then
       tmp = xlow
       xlow = xhigh
       xhigh = tmp
    endif

    do i = 1,mxiter

       xcur = (xlow + xhigh) / 2.0

       fxcur = f(xcur)

       if (fxcur <= 0d0) then
          xlow = xcur
       else
          xhigh = xcur
       endif

       diff = dabs(xhigh - xlow)

       if (diff <= toler) then
          bisect = xcur
          if (verbose) then
             write(6,"('iterations: ',i6,5f15.8)") &
                  i,xlow,xhigh,xcur,fxcur
          endif
          return
       endif
    enddo

   ! xcur = (xlow + xhigh)/2.0
   ! fxcur = f(xcur)

    write(6,"('Did not converge!!! iterations:  ',i6,5f15.8)") &
               niter,xlow,xhigh,xcur,fxcur
 

    bisect = xcur

  end function bisect

    double precision function newton(f, f_prime, x_init, mxiter, toler, verbose)

    real(ndp) :: x_init, xcur, xnew
    real(ndp), intent(in) :: toler
    integer(i4b), intent(in) :: mxiter
    logical :: verbose
    
    interface
       double precision function f(x)
         use constants
         real(ndp) :: x
       end function f
    end interface
    interface
       double precision function f_prime(x)
         use constants
         real(ndp) :: x
       end function f_prime
    end interface

  !  write(6,"('start ... bracket: ')")

    xcur = x_init
    
    do i = 1,mxiter

       xnew = xcur - f(xcur)/f_prime(xcur)

       fxcur = f(xnew)

       diff = dabs(xnew - xcur)
       xcur = xnew
       
       if (abs(fxcur) <= toler .or. diff <= toler ) then
          newton = xcur
          if (verbose) then
             write(6,"('iterations: ',i6,3f15.8)") &
                  i,xcur,fxcur,diff
          endif
          return
       endif

    enddo

    write(6,"('Did not converge!!! iterations:  ',i6,5f15.8)") niter,xcur,fxcur,diff

    newton = xcur

  end function newton

  !! Code provided by Tony Smith
   double precision function zbrent(fct,x1,x2,rtol,ftol,itmax,verbose)
! x1 and x2 must bracket the root.  rtol and ftol are convergence tolerances, the 
! first for the root itself and the second for the function value.
  !   implicit real(ndp) (a-h,o-z)
  !   implicit integer(i4b) (i-n)
     logical :: verbose
     
     integer(i4b) :: itmax
     integer(i4b) :: npos,nneg ! = 100
     real(ndp) :: x1, x2, rtol, ftol
     real(ndp) :: a, b, c, d,e, p, q, fa, fb, fc, xm, term2, tol1
     real(ndp), parameter :: toler = 1.0d-12
     real(ndp), parameter :: badval = -99.99d0
     real(ndp), parameter :: zero = 0d1     
     !      real(ndp), parameter :: eps = 3.0d-15

      interface
         double precision function fct(x)
           use constants
           real(ndp) :: x
         end function fct
      end interface
      
      tol1 = rtol

      npos = 0
      nneg = 0

      a = x1
      b = x2
      fa = fct(a)
      fb = fct(b)
   
      if (fb*fa > zero) then
         if ((fa > zero) .and. (fb > zero)) then
            npos = 1
            nneg = 0
          else
            npos = 0
            nneg = 1
         endif
         write(6,"(' Root must be bracketed in zbrent: ',4f15.8)") a,b,fa,fb
         zbrent = badval
         return
      endif
      c = b
      fc = fb
  
      do iter = 1,itmax
         if (fb*fc > zero) then
            c = a
            fc = fa
            d = b - a
            e = d
         endif
         if (dabs(fc) < dabs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
         endif
!         tol1 = two*eps*dabs(b) + 0.5d0*rtol
         xm = 0.5d0*(c-b)
         if ( (dabs(xm) <= tol1) .or. &
              (dabs(fb) <= ftol) ) then
            zbrent = b
            if (verbose) then
               write(6,"('iterations: ',i6,3f15.8)") &
                    iter,b,fb,xm
            endif
            return
         endif

         if ( (dabs(e) >= tol1) .and. &
              (dabs(fa) > dabs(fb)) ) then
            s = fb/fa
            if (dabs(c-a) <= toler) then
               p = two*xm*s
               q = one - s
            else
               q = fa/fc
               r = fb/fc
               p = s*(two*xm*q*(q-r)-(b-a)*(r-one))
               q = (q-one)*(r-one)*(s-one)
            endif
            if (p > zero) q = -q
            p = dabs(p)
            term1 = two*p
            term2 = dmin1(3.0d0*xm*q-dabs(tol1*q), &
                 dabs(e*q))
            ! & allows to continue code on the next line
            if (term1 < term2) then
               e = d
               d = p/q
            else
               d = xm
               e = d
            endif
         else
            d = xm
            e = d
         endif
         a = b
         fa = fb
         if (dabs(d) > tol1) then
            b = b + d
         else
            b = b + dsign(tol1,xm)
         endif
         fb = fct(b)
      enddo

      write(6,"(' zbrent exceeding maximum number of iterations')")

      zbrent = b

      return

   end function zbrent

! bla

  
end module rootfinding
