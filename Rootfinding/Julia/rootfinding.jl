function bisect(f, xlow, xhigh; mxiter::Int=30, toler::Real=1e-6, verbose=false)
    i = 1
    ## In order for the algorithm to work we need to bracket the root
    if f(xlow) * f(xhigh) > 0
        error("root not bracketed: f(xlow) = $(round(f(xlow),2)), f(xhigh) = $(round(f(xhigh),2))\n")
    ## We also want to have f(xlow) < 0 < f(xhigh) 
    elseif f(xlow) > 0
        tmp = xlow
        xlow = xhigh
        xhigh = tmp
    end
    
    for i = 1:mxiter
        #print("f(xlow) = $(round(f(xlow),2)), f(xhigh) = $(round(f(xhigh),2))\n")
        xcur = (xlow + xhigh)/2.0
        fxcur = f(xcur)
      
        if (fxcur <= 0.0)
            xlow = xcur
        else
            xhigh = xcur
        end
      
        diff = abs(xhigh - xlow)
        #@printf("Iterations: %4.i, %13.8f, %13.8f, %13.8f, %13.8f\n", i,  xlow, xhigh, xcur, fxcur)
        if diff <= toler
            if verbose
                @printf("Converged after %6.i iter, %13.8f, %13.8f, %13.8f\n",
                        i, xcur, fxcur, xhigh-xlow)
            end
            return xcur
        end
    end

   # xcur = (xlow + xhigh)/2.0
   # fxcur = f(xcur)

    error("Did not converge!  $i, val = $xcur,f(val) = $fxcur, Δ = $(xhigh-xlow)\n")   

end

function newton(f, fprime, x_init::Real,; mxiter::Int=30, toler::Real=1e-6, verbose=false)
    xcur = x_init
    i, fxcur, Δ = 1, 0.0, 0.0
    
    for i = 1:mxiter

        xnew = xcur - f(xcur)/fprime(xcur)

        fxcur = f(xnew)

        Δ = abs(xnew - xcur)
        xcur = xnew

       # @printf("iteration #: %6.i, %15.12f, %15.12f, %15.12f\n",  i,abs(xcur),fxcur,Δ)

        if (abs(fxcur) <= toler) & (Δ <= toler)
            if verbose
                @printf("Converged after %6.i iter, %13.8f, %13.8f, %13.8f\n",
                        i,xcur,fxcur,Δ)
            end
            return xcur
        end
    end
 
    error("Did not converge!  $i, val = $xcur,f(val) = $fxcur, Δ = $Δ\n")   
   
end


## Adapted from Fortran code provided by Tony Smith
function zbrent(fct, x1::Real, x2::Real;rtol::Real=1e-10,
                ftol::Real=1e-10,itmax::Int=40, verbose=false)
    ## x1 and x2 must bracket the root.  rtol and ftol are convergence tolerances, the 
    ## first for the root itself and the second for the function value.
    #     integer(i4b) :: npos,nneg
    toler = 1e-12
    eps = 3e-15
    d, e, xm = 0.0, 0.0, 0.0
    tol1 = rtol
    i = 1
    
    npos = 0
    nneg = 0

    a = x1
    b = x2
    fa = fct(a)
    fb = fct(b)
    if fb*fa > 0
        if (fa > 0) & (fb > 0)
            npos = 1
            nneg = 0
        else
            npos = 0
            nneg = 1
        end
        error("Root must be bracketed in zbrent: a = $a, b = $b, f(a) = $fa, f(b) = $fb")
    end
    c = b
    fc = fb

    for iter = 1:itmax
        if fb*fc > 0
            c = a
            fc = fa
            d = b - a
            e = d
        end
        if abs(fc) < abs(fb)
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        ##         tol1 = two*eps*dabs(b) + 0.5d0*rtol
        xm = 0.5*(c-b)
        if (abs(xm) <= tol1) | (abs(fb) <= ftol)
            if verbose
                @printf("Converged after %6.i iter, %13.8f, %13.8f, %13.8f\n",  i,b,fb,xm)
            end
            return b  
        end

        if (abs(e) >= tol1) | (abs(fa) > abs(fb))
            s = fb/fa
            if abs(c-a) <= toler
                p = 2.0*xm*s
                q = 1.0 - s
            else
                q = fa/fc
                r = fb/fc
                p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
                q = (q-1.0)*(r-1.0)*(s-1.0)
            end
            if p > 0
                q = -q
            end
            p = abs(p)
            term1 = 2.0*p
            term2 = min(3.0*xm*q-abs(tol1*q), abs(e*q))
            if term1 < term2
                e = d
                d = p/q
            else
                d = xm
                e = d
            end
        else
            d = xm
            e = d
        end
        a = b
        fa = fb
        if abs(d) > tol1
            b = b + d
        else
            b = b + copysign(tol1,xm) ####!!!! sign?
        end
        fb = fct(b)
    end

    error("zbrent exceeding maximum number of iterations $i, val = $b,f(val) = $fb, Δ = $xm\n")
end



### test functions

