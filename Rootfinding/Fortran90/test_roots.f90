module test_roots
  
  use constants
  use functions
  use rootfinding

  implicit none
contains
  
  subroutine test(time, res2, res3, res4, verbose, toler, N)
    
  use constants
  use functions
  use rootfinding

    
  real(ndp) :: num, x = 2.0
  real(ndp) :: xlow2, xhigh2, xlow3, xhigh3, toler
  integer(i4b) :: i, iterations, N
  logical :: verbose ! = .false.
  integer :: time_start, time_end, rate
  real(ndp) :: time
            
  real(ndp) :: res2(3), res3(2), res4(2)


  call system_clock(COUNT_RATE=rate)
  call system_clock(COUNT=time_start)
  
  do i = 1,N
     xlow2 = 0.001
     xhigh2 = 5
     iterations = 60
     !toler = 1d-10

     res2(1) = bisect(f2, xlow2, xhigh2, iterations, toler, verbose)
     res2(2) = newton(f2, f2p, xhigh2, iterations, toler, verbose)
     res2(3) = zbrent(f2, xlow2, xhigh2, toler, toler, iterations, verbose)
     if (verbose) then
        write(6,"('f2: 'f15.8)") sum(res2)/3
     endif

     xlow3 = 0.1
     xhigh3 = 1.99
     res3(1) = bisect(f3, xlow3, xhigh3, iterations, toler, verbose)
     res3(2) = zbrent(f3, xlow3, xhigh3, toler, toler, iterations, verbose)
     if (verbose) then
        write(6,"('f3: 'f15.8)") sum(res3)/2
     endif

     res4(1) = bisect(f4, xlow3, xhigh3, iterations, toler, verbose)
     res4(2) = zbrent(f4, xlow3, xhigh3, toler, toler, iterations, verbose)
     if (verbose) then
        write(6,"('f4: 'f15.8)") sum(res4)/2
     endif
  enddo
  
  call system_clock(COUNT=time_end)
  time = time_end - time_start
  time = time/rate
end subroutine test
  
end module test_roots
