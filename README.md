# Julia for Computational Economics
Provides sample code comparing Julia's performance with other languages used by economists.

[Julia](http://julialang.org) is a very young programming language that claims to be as easy as Matlab or Python while being as fast as C or Fortran.

In this repository I collect codes that I have written for graduate-level economics courses. All of these courses used programming languages other than Julia as their main language. I have used this opportunity to compare Julia's performance to other languages.

For now (spring 2017) I will post codes for my current course in Computational Methods for Economics, which teaches Fortran 90. Codes from previous courses will probably follow.

# Languages
- Julia
- Fortran
- (Matlab, Python, R, ...)

# Topics
- one dimensional root finding
- numerical differentiation
- (numerical optimization)
- (numerical integrations)
- (dynamic programming)

## Julia vs Fortran

Julia can be used like other intepreted languages (R, Python, Matlab). Type a command in a Julia console and you will get immediate response.

```julia
julia> 1 + 1
2
```

Fortran is a compiled language. That is, one has to write code to source file, compiles it, and then runs it.
For example, the file `demo-code.f90` contains the following code.
```fortran
program main
  integer :: x = 1+1              ! declare integer variable x to be 1+1
  write(6, "(i3)") x              ! print the x to the console (using 3 spaces)
end program main
```  
(Note that `!` starts a comment.) Then use a fortran compiler (e.g. `gfortran`) to compile the program.
```shell
gfortran -o demo.o demo-code.f90
```
This creates an executable file `demo.o` which can subsequently be run.
```
./demo.o                          # run code
  2
```
Luckily, it gives the same result as `Julia`.

### Speed comparison

When measuring speed, one does not want to measure compile time, but run time. For `Fortran` I will thus measure only the runtime of the executable. For `Julia` it is a bit more complicated. `Julia` uses "just in time compilation" (JIT). That is, when a function is called the first time it gets compiled automatically. Every subsequent time the function is called, `Julia` will automatically run the compiled version of the function. That is why I will call all functions once before I measure the runtime.

```julia
julia> f(x) = 1+1
f (generic function with 1 method)

julia> @time f(1)
  0.002551 seconds (330 allocations: 18.875 KB)
2

julia> @time f(2)
  0.000003 seconds (4 allocations: 160 bytes)
2
```
Furthermore I will measure the spead directly from `Julia`. One can send commands to the shell using ```run(`echo Hello`)```.
```julia
julia> @time run(`./a.out`)
  2
  0.009107 seconds (31 allocations: 1.313 KB)
```

