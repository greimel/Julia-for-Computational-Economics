
# Rootfinding

In this notebook I will do the speed comparisons.

We want to find the roots of the following real functions $f_i : \mathbb{R} \to \mathbb{R}$.

- $f_1(x) = f_1(x) = 0.5 x^{-0.5} + 0.5 x^{-0.2}$
- $f_2(x) = \log(a+x) - b$ for $(a, b) = 1.5, 1$
- $f_3(x) = - U'(a_0 + y_0 - x) + \beta (1+r) U'((1+r)\cdot x + y_1) $
for $U'(c) = c^{-\sigma}$ and $(a_0, y_0, y_1, r, \beta, \sigma) = (1, 1, 1, 0.03, 0.99, 2)$.

- $ f_4(x) = - U'(y_0 + a_0 - x) + \beta (1+r) \frac{1}{2} \bigl( U'(y_L + (1+r) x) + U'(y_H + (1+r) x))\bigr) $
where $U'(c) = c^{-\sigma}$, $(a_0, y_0, \bar y, r, \beta, \sigma) = (1, 1, 1, 0.03, 0.99, 2)$ and $\lambda \in \{ 0, 0.25, 0.5, 0.75 \}$.

These functions are found in `testfunctions.jl` and in the module `testfunctions.f90`.

We are using the methods 
1. bisection
2. newton
3. brent (Fortran Code by [Tony Smith](www.econ.yale.edu/smith/)).

The corresponding functions can be found in the files `rootfinding.jl` and `rootfinding.f90`.

## Description of the program
```
for i = 1:1000
    f_2 > bisect, newton, brent
    f_3 > bisect, brent
    f_4(λ) > bisect, brent λ = 0.25
end
```




```julia
include("Julia/testfunctions.jl")
include("Julia/rootfinding.jl")
```




    zbrent (generic function with 1 method)



The `Julia` version of the test program is provided below.


```julia
function test(; verbose=false)
    res2 = zeros(3)
    res3 = zeros(2)
    N = 1000
    
    ## Julia's task is a bit harder, the starting values vary randomly
    xlow2 = (rand(N) .- 0.5)./100 .+ 0.001
    xhigh2 = (rand(N) .- 0.5)./100 .+ 5.0
    
    xlow3 = (rand(N) .- 0.5)./100 .+ 0.01
    xhigh3 = (rand(N) .- 0.5)./100 .+ 1.99
    
    for i = 1:N 
        x = 2.0
        iterations = 40
        toler = 1e-10

        res2[1] = bisect(f2, xlow2[i], xhigh2[i], mxiter=iterations, toler=toler,verbose=verbose)
        res2[2] = newton(f2, f2p, xhigh2[i], mxiter=iterations, toler=toler, verbose=verbose)
        res2[3] = zbrent(f2, xlow2[i], xhigh2[i], rtol=toler, ftol=toler, itmax=iterations, verbose=verbose)
                
        res3[1] = bisect(f3, xlow3[i], xhigh3[i], mxiter=iterations, toler=toler, verbose=verbose)
        res3[2] = zbrent(f3, xlow3[i], xhigh3[i], rtol=toler, ftol=toler, itmax=iterations, verbose=verbose)
        
        res4 = zeros(2)

        res4[1] = bisect(x::Real -> f4(x, 0.25), xlow3[i], xhigh3[i], mxiter=iterations, toler=toler, verbose=verbose)
        res4[2] = zbrent(x::Real -> f4(x, 0.25), xlow3[i], xhigh3[i], rtol=toler, ftol=toler, itmax=iterations, verbose=verbose)
        
        if verbose
            mean2 = mean(res2)
            mean3 = mean(res3)
            mean4 = mean(res4)
            print("$mean2, $mean3, $mean4, $xlow2, $xhigh2\n")
        end
    end
end

# warm up
test()
```

    WARNING: Method definition test() in module Main at In[9]:2 overwritten at In[11]:2.
    WARNING: Method definition #test(Array{Any, 1}, Main.#test) in module Main overwritten.


I also provide the main program of the `Fortran` version (see also `test.f90`) for convenience.

```fortran
program main
  
  use constants
  use functions
  use rootfinding
 
  real(ndp) :: x = 2.0
  real(ndp) :: xlow2, xhigh2, xlow3, xhigh3, toler
  integer(i4b) :: i, iterations
  logical :: verbose = .false.
  
  real(ndp) :: res2(3), res3(2), res4(2)

  do i = 1,1000
     xlow2 = 0.001
     xhigh2 = 5
     iterations = 60
     toler = 1d-10

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

  end program main
```

Call the `makefile` to compile the code.


```julia
; cd Fortran90
```

    /Users/Fabian/Dropbox/Yale_Courses/Comp_Econ/Julia-for-Computational-Economics/Rootfinding/Fortran90



```julia
; make
```

    make: `test.out' is up to date.



```julia
; cd ..
```

    /Users/Fabian/Dropbox/Yale_Courses/Comp_Econ/Julia-for-Computational-Economics/Rootfinding



```julia
M = 100
time_julia = zeros(M)
time_fortran = zeros(M)

for i = 1:M
    time_fortran[i] = @elapsed run(`./Fortran90/test.out`)
    time_julia[i] = @elapsed test()
end
```


```julia
using Plots, PyPlot
```


```julia
pyplot()

histogram(time_fortran[2:end], label="Fortran", alpha=0.5)
histogram!(time_julia[2:end], label="Julia", alpha=0.5)

plot!([median(time_fortran);median(time_julia)], linetype=:vline, label="medians")
```




<img src="data:image/png;base64,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IdAsS0pKNGvWLEVHR2v8+PHasWOHY8MCAICux+fz6f7779f48eM1bdo0FRUV6YMPPnB7rGZdsoCuWLFCoaGhOnz4sCTp5MmTkqRf/OIXiouL05tvvqn3339fd9xxhwoKChQWFubMxAAAoMPq1q2bampqAvcrKipafM4jjzyiIUOGaO/evQoJCdH3v/99S89zU5OfYPl8Pr300kv6zW9+E3hs6NChkqRNmzZpyZIlkqTJkydrxIgRevfddx0YtWvKzs52ewQAFpBVoHl+v1+S9I1vfEPvvfeeJOno0aOWsnP69GldccUVCgkJ0UcffaS33347qLPaoclPsPLz8zVw4ED95je/0TvvvKOIiAilpKQoJiZG1dXVgbIlSZGRkSoqKnJs4K5m1apV+va3v+32GABaQFZhutDL+tt63FToZf1btX11dbXCw8OVnJysuXPnasKECRo7dqymTZsW2ObiMVcN/eu//qvmz5+v3/3udxo1apS++93vtmt2JzT5s66urlZhYaGuvfZa/du//Zs++OADxcfH68CBA4EGCmds2LDB7REAWEBWYbrWnPFnt+LiYpWVlcnj8Sg8PFy7du1qcruXXnop8OeUlJTAn2NiYrR///6gz2mnJr8i9Hg8CgsL07x58yR9tWMjR47U/v371b1798DxWJJUUFAgj8fT7JvMnj1bXq+33i0uLq7R6vOff/653sj6Y5OvkffRR/XuHy8u1htZf1RFZWW9x1NSUpSenl7vsaKiInm9XuXl5dV7/Nlnn1VSUlK9x3w+n7xeb6OPLDMzM5u8zs3cuXMb7ce2bdsCp5/WtXTpUq1fv77eYzk5OfJ6vSotrX9mx8X9iIiI6BT7URf7wX50xv24mNWOvh8XsR8ddz+ef/75Rs9309NPP60bb7xRTz31lMLDw90exzZ33313oNMsX7680d+H+C/xkdTMmTP10EMPadasWTp69Kiuu+467d27V//yL/+iyMhIpaSkaPfu3ZozZ84lD3LPycnRpEmTtGfPHsXGxjY7aHFxsdauTNP8aWM1bODXHzvmV5/Ww1/+RU/3/a5GdRtQ7znHT53Ry+8d1NJHUjR8+HBLPxAAADqz1vy/F63X1M+3qccu+XXsc889p/vvv18rVqxQWFiYXnjhBQ0fPlyPP/645s+fr+joaPXs2VOvvPIKZxACAADUccmCdfGiXg0NHTpUW7duDepQ+FpSUpKeeOIJt8cA0AKyCqAu1iI0XEvHtwEwA1kFUBcFy3DLli1zewQAFpBVAHVRsAAAAGxGwQIAAB3OwYMHFRUVJemrKxF85zvfcXmi+ihYhmt47RUAZiKrgPMuXvl9+PDhxi3bR8EyXHJystsjALCArAKXFhoaqpUrV2ratGm66qqr9MYbb+jxxx/XlClTNHr0aP3tb38LbLtt2zZdf/31mjJliqZNm6b//d//DfxdamqqoqOjNWXKlHqrJxQWFmrAgK+vlXnPPfdo6tSpiomJ0a233hq4QPrF7VJTUzV58mRFR0frrbfekvTVotN33323rr32Wk2cOFEzZ85s1z7buSwRgmDNmjVujwDAArKKjuTzsuMqv3CuTc/t0723Lr9sWKuf17dvX7333nvavn27brvtNv3nf/6ndu/erddee00///nPtWvXLh09elSpqanatm2b+vTpo/z8fF1//fUqLCzUtm3btHnzZuXm5qp3796aP39+vdevu47hM888o0GDBkmS0tPTlZKSoueee06SdPbsWcXExCg1NVVbt27VQw89pLy8PL311ls6e/asDhw4IEk6c+ZMm34+F1GwDMep30DHQFbRUZyp+FI//NMDqvXXtun5YSGhen3O79Q/vG+rnnfXXXdJkiZPniyfz6e5c+dKkqZOnapPPvlEkvTWW28pPz9fN9xwQ2Dt427duqmoqEjbt2/XXXfdpd69e0uSFi9erJ07dzb5XhkZGcrIyFBFRYUqKys1ePDgwN/16tVLt99+uyQpLi5OR44ckSRNmDBBhw4dUmJiom644QbNnj27VfvXEAULAIAupH94X71y63Pt+gSrteUqJCQksA7hxdVfevToEbhfXV0tSfL7/br55puVkZHRptkkKTs7W88++6z+8Y9/aNCgQfrTn/5Ub+Honj17Bv4cFhammpoaSV9dYP3DDz/U9u3b9fbbbys5OVl79+5Vv3792jQHBQsAgC6mLV/xtUfDZY8vdX/GjBn69a9/rf3792vcuHGSpN27d2vKlCm66aabtGLFCj388MPq3bu3XnzxxSZf48yZM+rbt68GDBigqqoqrVu3ztJ7f/bZZxowYIC+973vacaMGXrjjTd07NixNhcsDnI3XMOV3QGYiawCl1b3+Kjm7o8aNUp/+MMftHjxYk2cOFFjx47VM888I0maNWuW7rzzTsXGxmrq1KmKjIxs8jVmzpyp6OhojR49Wt/5znc0ceJES++9f/9+fetb39LEiRMVGxure++9V9dee22b95lPsAzn8/ncHgGABWQVuLSLX8NJUu/evevdHzFihL788svA/RtvvFH/93//1+TrPProo3r00UcD93/9619LkiIjI3Xq1ClJXx2zVfcMQ0l67LHHGm3XcJaZM2e2+8zBuvgEy3BpaWlujwDAArIKoC4KFgAAgM0oWAAAADajYBmutLTU7REAWEBWAdTFQe6GW7RokbZs2eL2GABaQFZhmkOHDrk9Qqdk9edKwTJcamqq2yMAsICswhSDBw9WRESE7rnnHrdH6bQiIiLqXR2+KRQsw8XGxro9AgALyCpM4fF4dOjQoWa/ti48e0yP7VytX33rp4rsd6WD03UOgwcPbnF5LAoWAACdjMfjabYA9DnVT5d93l9jJ1yr6IGjHJys6+AgdwAAAJtRsAy3fv16t0cAYAFZBcznZE4pWIbLyclxewQAFpBVwHxO5pSCZbi1a9e6PQIAC8gqYD4nc0rBAgAAsBkFCwAAwGYULAAAAJtRsAzn9XrdHgGABWQVMJ+TOaVgGS4xMdHtEQBYQFYB8zmZUwqW4eLj490eAYAFZBUwn5M5pWABAADYjIIFAABgMwqW4bKystweAYAFZBUwn5M5pWAZLjMz0+0RAFhAVgHzOZlTCpbhNm7c6PYIACwgq4D5nMwpBQsAAMBmFCwAAACbUbAAAABsRsEy3MKFC90eAYAFZBUwn5M5pWAZjqtDAx0DWQXMx5XcEZCQkOD2CAAsIKuA+ZzMKQULAADAZhQsAAAAm1GwDJedne32CAAsIKuA+ZzMKQXLcKtWrXJ7BAAWkFXAfE7mlIJluA0bNrg9AgALyCpgPidzSsEyXEREhNsjALCArALmczKnFCwAAACbUbAAAABsRsEyXFJSktsjALCArALmczKnFCzDeTwet0cAYAFZBcznZE4pWIZbtmyZ2yMAsICsAuZzMqcULAAAAJtRsAAAAGxGwTJcXl6e2yMAsICsAuZzMqcULMMlJye7PQIAC8gqYD4nc0rBMtyaNWvcHgGABWQVMJ+TOaVgGY5Tv4GOgawC5uMyDQAAAB0YBQsAAMBmFCzDpaenuz0CAAvIKmA+J3NKwTKcz+dzewQAFpBVwHxO5pSCZbi0tDS3RwBgAVkFzOdkTilYAAAANqNgAQAA2IyCZbjS0lK3RwBgAVkFzOdkTilYhlu0aJHbIwCwgKwC5nMypxQsw6Wmpro9AgALyCpgPidzSsEyXGxsrNsjALCArALmczKnFCwAAACbUbAAAABsRsEy3Pr1690eAYAFZBUwn5M5pWAZLicnx+0RAFhAVgHzOZlTCpbh1q5d6/YIACwgq4D5nMwpBQsAAMBmFCwAAACbUbAAAABsRsEynNfrdXsEABaQVcB8TuaUgmW4xMREt0cAYAFZBcznZE4pWIaLj493ewQAFpBVwHxO5pSCBQAAYDMKFgAAgM0oWIbLyspyewQAFpBVwHxO5pSCZbjMzEy3RwBgAVkFzOdkTilYhtu4caPbIwCwgKwC5nMypxQsAAAAm1GwAAAAbEbBAgAAsBkFy3ALFy50ewQAFpBVwHxO5pSCZTiuDg10DGQVMB9XckdAQkKC2yMAsICsAuZzMqcULAAAAJtRsAAAAGxGwTJcdna22yMAsICsAuZzMqcULMOtWrXK7REAWEBWAfM5mVMKluE2bNjg9ggALCCrgPmczCkFy3ARERFujwDAArIKmM/JnFKwAAAAbEbBAgAAsBkFy3BJSUlujwDAArIKmM/JnFKwDOfxeNweAYAFZBUwn5M5pWAZbtmyZW6PAMACsgqYz8mcUrAAAABsRsECAACwGQXLcHl5eW6PAMACsgqYz8mcNluwfvvb3yo0NFRbtmyRJJWUlGjWrFmKjo7W+PHjtWPHDkeG7MqSk5PdHgGABWQVMJ+TOb1kwSosLNR//dd/KS4uLvDYL37xC8XFxenw4cN66aWXNG/ePNXU1DgyaFe1Zs0at0cAYAFZBcznZE6bLFh+v18/+tGPtGbNGvXo0SPw+KZNm7RkyRJJ0uTJkzVixAi9++67zkzaRXHqN9AxkFXAfK5fpmH16tW6/vrrNXHixMBjp06dUnV1tYYOHRp4LDIyUkVFRcGfEgAAoAPp1vCBgwcPavPmzRxfBQAA0EaNPsHasWOHCgsLdfXVVysqKkrvvfeefvKTn2jTpk3q1q2bTp48Gdi2oKDA0sdts2fPltfrrXeLi4tTVlZWve0+//xzvZH1xyZfI++jj+rdP15crDey/qiKysp6j6ekpCg9Pb3eY0VFRfJ6vY3OHnj22WcbXTbf5/PJ6/UqOzu73uOZmZlauHBho7nmzp3baD+2bdsmr9fbaNulS5dq/fr19R7LycmR1+tVaWlpk/tRd1868n7UxX6wH51xPy7O09H34yL2o2vsx4YNGzrFflj975Gent7u/cjMzFRUVJRiYmICnWb58uWNXi/E7/f7Gz1ax/Tp0/XTn/5Ut956qxYtWqTIyEilpKRo9+7dmjNnjgoKChQWFtbkc3NycjRp0iTt2bNHsbGxzb2NiouLtXZlmuZPG6thA/sHHs+vPq2Hv/yLnu77XY3qNqDec46fOqOX3zuopY+kaPjw4c2+fkeVkpKitLQ0t8cA0AKyio7k8Kl8/fjNn+rFWasVPXCU2+M4Jlg5barvNPqKsKGQkBBd7GCPP/645s+fr+joaPXs2VOvvPLKJcsV7ME/2EDHQFYB8zmZ0xYL1vbt2wN/Hjp0qLZu3RrUgQAAADo6ruQOAABgMwqW4Roe9AfATGQVMJ+TOaVgGW7RokVujwDAArIKmM/JnFKwDJeamur2CAAsIKuA+ZzMKQXLcC1d3gKAGcgqYD4nc0rBAgAAsBkFCwAAwGYULMM1XFIAgJnIKmA+J3NKwTJcTk6O2yMAsICsAuZzMqcULMOtXbvW7REAWEBWAfM5mVMKFgAAgM0oWAAAADajYAEAANiMgmU4r9fr9ggALCCrgPmczCkFy3CJiYlujwDAArIKmM/JnFKwDBcfH+/2CAAsIKuA+ZzMKQULAADAZhQsAAAAm1GwDJeVleX2CAAsIKuA+ZzMKQXLcJmZmW6PAMACsgqYz8mcUrAMt3HjRrdHAGABWQXM52ROKVgAAAA2o2ABAADYjIIFAABgMwqW4RYuXOj2CAAsIKuA+ZzMKQXLcFwdGugYyCpgPq7kjoCEhAS3RwBgAVkFzOdkTilYAAAANqNgAQAA2IyCZbjs7Gy3RwBgAVkFzOdkTilYhlu1apXbIwCwgKwC5nMypxQsw23YsMHtEQBYQFYB8zmZUwqW4SIiItweAYAFZBUwn5M5pWABAADYjIIFAABgMwqW4ZKSktweAYAFZBUwn5M5pWAZzuPxuD0CAAvIKmA+J3NKwTLcsmXL3B4BgAVkFTCfkzmlYAEAANiMggUAAGAzCpbh8vLy3B4BgAVkFTCfkzmlYBkuOTnZ7REAWEBWAfM5mVMKluHWrFnj9ggALCCrgPmczCkFy3Cc+g10DGQVMB+XaQAAAOjAKFgAAAA2o2AZLj093e0RAFhAVgHzOZlTCpbhfD6f2yMAsICsAuZzMqcULMOlpaW5PQIAC8gqYD4nc0rBAgAAsBkFCwAAwGYULMOVlpa6PQIAC8gqYD4nc0rBMtyiRYvcHgGABWQVMJ+TOaVgGS41NdXtEQBYQFYB8zmZUwqW4WJjY90eAYAFZBUwn5M5pWABAADYjIIFAABgMwqW4davX+/2CAAsIKuA+ZzMKQXLcDk5OW6PAMACsgqYz8mcUrAMt3btWrdHAGABWQXM52ROKVgAAAA2o2ABAADYjIIFAABgMwqW4bxer9sjALCArALmczKnFCzDJSYmuj0CAAvIKmA+J3NKwTJcfHy82yMAsICsAuZzMqcULAAAAJtRsAAAAGxGwTJcVlaW2yMAsICsAuZzMqcULMNlZma6PQIAC8gqYD4nc0rBMtzGjRvdHgGABWQVMJ+TOaVgAQAA2IyCBQAAYDMKFgAAgM0oWIZbuHCh2yMAsICsAuZzMqcULMNxdWigYyCrgPm4kjsCEhIS3B4BgAVkFTCfkzmlYAEAANiMggUAAGAzCpbhsrOz3R4BgAVkFTCfkzmlYBlu1apVbo8AwAKyCpjPyZxSsAy3YcMGt0cAYAFZBcznZE4pWIaLiIhwewQAFpBVwHxO5pSCBQAAYDMKFgAAgM0oWIZLSkpyewQAFpBVwHxO5pSCZTiPx+P2CAAsIKuA+ZzMKQXLcMuWLXN7BAAWkFXAfE7mlIIFAABgMwoWAACAzShYhsvLy3N7BAAWkFXAfE7mlIJluOTkZLdHAGABWQXM52ROKViGW7NmjdsjALCArALmczKnFCzDceo30DGQVcB8XKYBAACgA6NgAQAA2IyCZbj09HS3RwBgAVkFzOdkTilYhvP5fG6PAMACsgqYz8mcUrAMl5aW5vYIACwgq4D5nMwpBQsAAMBmFCwAAACbUbAMV1pa6vYIACwgq4D5nMwpBctwixYtcnsEABaQVcB8TuaUgmW41NRUt0cAYAFZBcznZE4pWIaLjY11ewQAFpBVwHxO5rTJglVZWak77rhD11xzjSZOnKgZM2YoPz9fklRSUqJZs2YpOjpa48eP144dOxwbFgAAoCO45CdYixcvVl5ennJzc+X1evWjH/1IkrRixQrFxcXp8OHDeumllzRv3jzV1NQ4NjAAAIDpmixYPXv21MyZMwP3p02bpsLCQknSq6++qiVLlkiSJk+erBEjRujdd991YNSuaf369W6PAMACsgqYz8mcWjoG65lnntHtt9+uU6dOqbq6WkOHDg38XWRkpIqKioI2YFeXk5Pj9ggALCCrgPmczGm3ljZYuXKl8vPz9cILL7DWlgvWrl3r9ggALCCrgPmczGmzn2A9+eSTysrK0ltvvaXw8HANHDhQ3bp108mTJwPbFBQUyOPxNPsms2fPltfrrXeLi4tTVlZWve0+//xzvZH1xyZfI++jj+rdP15crDey/qiKysp6j6ekpDRaLbuoqEher1d5eXn1Hn/22WeVlJRU7zGfzyev16vs7Ox6j2dmZmrhwoWN5po7d26j/di2bZu8Xm+jbZcuXdro48mcnBx5vd5GFz9jP9gP9oP9YD/Yj2Dvx4YNGzrFfjj53yMzM1NRUVGKiYkJdJrly5c3er0Qv9/vb/SopNWrV+sPf/iD/vKXv6hfv36BxxctWqTIyEilpKRo9+7dmjNnjgoKChQWFtboNXJycjRp0iTt2bOnxVMji4uLtXZlmuZPG6thA/sHHs+vPq2Hv/yLnu77XY3qNqDec46fOqOX3zuopY+kaPjw4c2+PgAA+MrhU/n68Zs/1YuzVit64Ci3x+nwmuo7TX5F+Nlnn+nnP/+5Ro0apenTp8vv9ys8PFx///vf9fjjj2v+/PmKjo5Wz5499corrzRZrgAAALqqJgvWiBEjVFtb2+QThg4dqq1btwZ1KHzN6/Vqy5Ytbo8BoAVkFTCfkznlSu6GS0xMdHsEABaQVcB8TuaUgmW4+Ph4t0cAYAFZBcznZE4pWAAAADajYAEAANiMgmW4htfkAGAmsgqYz8mcUrAMl5mZ6fYIACwgq4D5nMwpBctwGzdudHsEABaQVcB8TuaUggUAAGAzChYAAIDNKFgAAAA2o2AZrqkVvwGYh6wC5nMypxQsw3F1aKBjIKuA+biSOwISEhLcHgGABWQVMJ+TOaVgAQAA2IyCBQAAYDMKluGys7PdHgGABWQVMJ+TOaVgGW7VqlVujwDAArIKmM/JnFKwDLdhwwa3RwBgAVkFzOdkTilYhouIiHB7BAAWkFXAfE7mlIIFAABgMwoWAACAzShYhktKSnJ7BAAWkFXAfE7mlIJlOI/H4/YIACwgq4D5nMwpBctwy5Ytc3sEABaQVcB8TuaUggUAAGAzChYAAIDNKFiGy8vLc3sEABaQVcB8TuaUgmW45ORkt0cAYAFZBcznZE4pWIZbs2aN2yMAsICsAuZzMqcULMNx6jfQMZBVwHxO5rSbY+9kkLKyMpWXl7fqOREREerXr1+QJgIAAJ1JlytYZWVleuKxNNWeb13B6tm3vx5M/iUlCwAAtKjLFazy8nLVni/XzdEjNGxgf0vPKT37pd48dEw+n8/xgpWenq4VK1Y4+p6o36HuAAAOR0lEQVQAWo+sAuZzMqddrmBdNGxgf8sFy00+n8/tEQBYQFYB8zmZUw5yN1xaWprbIwCwgKwC5nMyp132EywAALqqc+fOSZJKS0t1WWUEJ3IFAQULAIAupKysTM+99LL0DenZ37+mnucjNLh3d/3qZ8soWTbiK0LDlZaWuj0CAAvIKjqK8vJylV346s/94u5Qrym3q/TchS5xHKGTOaVgGW7RokVujwDAArKKjqh3/8GKGDDE7TEc42ROKViGS01NdXsEABaQVcB8TuaUgmW42NhYt0cAYAFZBcznZE4pWAAAADajYAEAANiMgmW49evXuz0CAAvIKmA+J3NKwTJcTk6O2yMAsICsAuZzMqcULMOtXbvW7REAWEBWAfM5mVMKFgAAgM1YKgcAgC7uQlWlTpw4EbjP2oTtR8ECAKALq/KVKzc3VyvX1SoiopcksTahDfiK0HBer9ftEQBYQFbRUVWdP6eq0HD1mHiLBsX/uFOvTehkTvkEy3CJiYlujwDAArKKjq73oGG6bPBwSdJ5l2cJFidzyidYhouPj3d7BAAWkFXAfE7mlIIFAABgMwoWAACAzShYhsvKynJ7BAAWkFXAfE7mlIJluMzMTLdHAGABWQXM52ROKViG27hxo9sjALCArALmczKnFCwAAACbcR0sAAAMUFZWpvLy8sB9lqvp2ChYAAC4rKysTL9a+aTOVH39GMvVdGx8RWi4hQsXuj0CAAvIKtqjvLxcZ6qknrHf6/TL1bjJyZzyCZbhuDo00DGQVdihz+DhnX65GjdxJXcEJCQkuD0CAAvIKmA+J3NKwQIAALBZh/+KsLKqSidOnLC8/YkTJ1R14UIQJwIAAF1dhy5Y5ecrlJubq9o1T6tXr16WnlN27pyOHM5T1dQxQZ7OHtnZ2fr2t7/t9hgAWkBWAfM5mdMOXbB8FZXqqVrdfPXlirrickvPySv6TJ8c2KcL1dVBns4eq1at4h9toAMgq4D5nMxphy5YFw0d0E/DBva3tO3J02eDPI29NmzY4PYIACwgq4D5nMwpB7kbLiIiwu0RAFhAVgHzOZlTChYAAIDNOsVXhAAAmM7NtQbrvveJEyd0oYWz6S9UVdY7Q591EVuPgmW4pKQkPfHEE26PAaAFZBXNcXOtwYbv7Ssv0+FTRzVgfNPbV/nKlZubq5XrahUR0cvRWYPNyZxSsAzn8XjcHgGABWQVzam71mCfwcPlO12i0t1Z8vl8QS8tDd+79siHuvD22ktuX3X+nKpCw9Vj4i0aFDnK0VmDzcmcUrAMt2zZMrdHAGABWYUVbq41ePG9y7+wdnHu3oOGdbp1EZ3MKQe5AwAA2IxPsAAAcEHdA8lPnDihqqqqFp6BjoSCZbi8vDxdc801bo8BoAVkFa3R8EByX3mZDn58RANvqtRlbg/XiTmZU74iNFxycrLbIwCwgKyiNeodSB7/Y4XHzFRltb/DLOPWUTmZUwqW4dasWeP2CAAsIKtoi4sHkvfqN8jtUboEJ3NKwTIcp34DHQNZBcznZE4pWAAAADbjIHcAADqgusvfsJSNefgEy3Dp6elujwDAArIKJ11c/uaXTz6vXz75vB576lmdPXvW7bGM52ROKViG8/l8bo8AwAKyCifVXf6m15TbVXruAr+DFjj5M+IrQsOlpaW5PQIAC8gq3NCnky1lE2xO5pRPsAAAAGxGwQIAALAZXxEarrS0VIMHD3Z7DAAtIKswSd0zDIO9zmHd95LMPqPRyZxSsAy3aNEibdmyxe0xALSArMIUF88wPPP/O1Uw1zls+F6SNLh3d/3qZ8uMLFlO5pSCZbjU1FS3RwBgAVmFKeqeYdhn8HDVHvlQlYfWBmWdw4bv5TtdotLdWfL5fEYWLCdzSsEyXGxsrNsjALCArMI0fQYP12WDh6v8ixOOvZdk9hmNTuaUg9wBAABsxidYACxpeCBrS0w+0BXoCC5UVerEia8/faqpqVFYWJikxgeu1902GAe1N5ylPfnuSAfFtwcFy3Dr16/X/fff7/YY6OLKysr0xGNpqj1vvWD17NtfDyb/slP+w9kUsgo7VfnKlZubq5XrahUR0UsXqip1aP9eXTMuRj169Kh34HpI7YV629p9UHvDWaSvD2RvLbcPincypxQsw+Xk5PCPNlxXXl6u2vPlujl6hIYN7N/i9qVnv9Sbh44Ze6BrMJBV2Knq/DlVhYarx8RbNChylEqOfKgv9x5S6LiZGhQ5qv6B65X1t7X7oPaGs9Q9kL213D4o3smcUrAMt3btWrdHAAKGDexvqWB1RWQVwdB70LB6B6o3vN/ctsGaRWr/gexuHRTvZE45yB0AAMBmFCwAAACb8RWhRZVVVfXOoLCi7hkfVnXWsykAAF1Da884bO32HeUsRAqWBeXnK5Sbm6vaNU+rV69elp5TWVWlfQc/1IRrx6pH9+6W36vhmVder5flN4AOgKwCzZ9x2FQJau327T0L0cmcUrAs8FVUqqdqdfPVlyvqisstPSev6DN9mLtHN0YNtfycps68SkxMbPPcAJxDVoHmzzhssmC1cvv2noXoZE4pWK0wdEA/y2dQnTx9ttXPaUp8fHybnwvAOWQV+Fprzzhs7fZtPQvRyZxykDsAAIDNKFgAAAA24ytCw2VlZen22293ewwYrLVrBErOnHXTljNvnTobKBg/M7LamKm/m2i/umf+2b3uYUN1f4+aWmex7iwNz95v+PvU3py25gxGCpbh0tPT+Ucbl9SWNQKl4K8T2JYzb52YSwrez4ys1mfq7ybar+6Zf6qttnXdw4YanjXYcJ3FurN07xZab71GqfEZhu3JaXNnMDalTQXrk08+0YIFC1RaWqr+/fvrv//7vzVmzJg2DYzmDRkyxO0RYLDWrhEoObNOYFvOvHVq/cJg/czIan2m/m6i/eqe+RdaU2nruocNNTxrsOE6iw1nqbteY1NnGLYnp82dwdiUNhWsxYsXa8mSJZo/f742b96sBQsWaNeuXW0eGkD7mLpGYHvPog0mU39mnQ0/586r96BhUmXrPqFsq4tnDV5qncW6s9i5ZmJzs7T0+q0+yL2kpER79uzRD3/4Q0nS97//fR07dkxHjhxp06AAAACdTasL1rFjxzR8+HCFhn79VI/Ho6KiIlsHAwAA6KiCepD7+fNffXh26NChFrctKSlR8ckS7dyfp/59egceP9mjQhoh7Tr0iQqqwus9p+D4SZ0959N7H36sgpOnLM1k8nPOlJ9TfmGRNm/erAEDBkiSdu7cqVdeeeWSzwkNDVVtba2l1+c5bX+OqXOdPn1aBcc+1c7eYfVy05ymfs/sfh+7fv+tMOVn1lRW+d1s/c+5+GSJ9u3bp+Li4lbNZ7qSkhJ9caJY5Tl/U/hlA3T6syOqLD+jT/fu1LniI+26rwvnW/3cKl+ZpMt0bN97OldUFrT3rig7rbMFnwSyc/r0aX1WWHDJn0Nrt29uloavJTXOacPf6+buN5ylouy0Kk8Ua9++fSotLZX0de+RJPlb6eTJk/5+/fr5a2pqAo8NGzbMn5+f32jbjIwMvyRu3Lhx48aNG7dOf8vIyAh0oFZ/gjVkyBDFxsbq5Zdf1oIFC/Taa6/pyiuv1FVXXdVo2xkzZigjI0MjR45s1anaAAAAHcX58+dVUFCgGTNmBB4L8fv9/ta+0OHDh3Xffffpiy++UL9+/fTb3/5WY8eOtXVYAACAjqpNBQsAAACXxlqEAAAANqNgueCTTz7Rt771LY0ePVrXXXfdJc+y/POf/6wxY8Zo9OjRuvPOOwPrHxUXF2vmzJkaM2aMYmJi9IMf/EBffPGFk7sAdAntzWpdKSkpCg0N1b59+4I9NtDl2JHVM2fO6J577tHo0aM1btw4PfLII+0bqrVnEaL9brzxRv/vf/97v9/v97/22mv+KVOmNNqmvLzc/0//9E/+w4cP+/1+vz8xMdGflJTk9/v9/hMnTvh37twZ2DYpKcl/3333OTA50LW0N6sX7dq1yz979mx/VFSUf+/evcEfHOhi7MjqHXfc4V+9enXg/okTJ9o1EwXLYVYvc/Hqq6/6Z82aFbj/4Ycf+q+44oomX/O1117zT58+PTgDA12UXVn1+Xz+qVOn+j/99FP/yJEjKViAzezI6scff+z3eDy2zsVXhA6zeiX8oqIiRUZGBu6PHDlSx48fb3Shv9raWq1Zs6bNq4MDaJpdWU1OTtbSpUs1YsQIZwYHuhg7snro0CGNGDFCS5Ys0eTJkzVz5kx98MEH7ZqLgtXBPfDAAxo4cKAefPBBt0cB0MDbb7+twsJC3XvvvW6PAqAZ1dXV2rVrl+bNm6f3339fy5cv1/e+9z3V1NS0+TUpWA678sorVVxcXO+TqKKiInk8nnrbeTweFRQUBO4fPXq0UUN/8MEH9fnnn2vTpk1BnxvoauzI6l//+lfl5ubqqquuUlRUlD799FPNnj1b//M//+PUbgCdnh1Z9Xg8uuKKK3TDDTdIkmbOnKmqqioVFha2eS4KlsPqXglf0iWvhD9z5kzl5ubq8OHDkqTnnntOd999d+DvH3zwQeXn5+v1119XWFiYczsAdBF2ZHXlypU6duyYjhw5oqNHj+qKK67Qm2++qVtuucXZnQE6MTuyOmnSJPXt21f79++XJO3atUvSV+WtzWw9oguWfPTRR/64uDh/dHS0f8qUKf6DBw/6/X6//9FHH/WvW7cusN2f/vQn/zXXXOO/+uqr/XfccYf/yy+/9Pv9fv/OnTv9oaGh/m9+85v+mJgYf0xMjH/OnDmu7AvQmbU3qw1xFiEQHHZkNScnx3/dddf5J0yY4J86dap/x44d7ZqJK7kDAADYjK8IAQAAbEbBAgAAsBkFCwAAwGYULAAAAJtRsAAAAGz2/wAXk2Jo+ioqhwAAAABJRU5ErkJggg==" />


