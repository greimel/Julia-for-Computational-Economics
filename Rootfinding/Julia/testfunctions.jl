f1(x::Real; a::Real=0.5, b::Real=0.5) = a*x^(-0.5) + b*x^(-0.2)
f1p(x::Real; a::Real=0.5, b::Real=0.5) =
    - 0.5*a*x^(-1.5) - 0.2 * b*x^(-1.2)

f2(x::Real; a::Real=1.5, b::Real=1.0) =  log(a + x) - b
f2p(x::Real; a::Real=1.5, b::Real=1.0) = 1.0/(a+x)

function f3(x::Real; a0=1, y0=1, y1=1, r=0.03, β=0.99, σ=2)
    Up(c) = c^(-σ)
    -Up(a0 + y0 - x) + β * (1+r) * Up((1+r)*x + y1)
end

function f4(x::Real, λ::Real; a0=1, y0=1, ybar=1, r=0.03, β=0.99, σ=2)
    Up(c) = c^(-σ)
    yL = (1-λ) * ybar
    yH = (1+λ) * ybar
    -Up(a0 + y0 - x) + β*(1+r)/2 * (Up((1+r)*x + yH) + Up((1+r)*x + yL))
end
