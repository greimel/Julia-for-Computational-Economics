## Julia has a fantastic `Parameters.jl` package. It makes the code cleaner, without losing performance.

using Parameters

type Par{T <: Real}
    α::T
    β::T
end

### Implementation 1: unpack in outer function
function f1(x::T, y::T, α::T, β::T) where T <: Real
    x^α * y^β
end

function g1(xx::Array{T}, yy::Array{T}, α::T, β::T) where T <: Real
    mean(f1.(xx, yy, α, β))
end

function h1(n::Int, p::Par{Float64})
    @unpack α, β = p
    out = 0.0
    for i in 1:n
        xx = rand(500, 500)
        yy = rand(500, 500)
        out += g1(xx, yy, α, β)/n
    end
end

### Implementation 2: unpack in inner function
## inlining makes a difference!
@inline function f2(x::T, y::T, p::Par{T}) where T <: Real
    @unpack α, β = p
    x^α * y^β
end

function g2(xx::Array{T}, yy::Array{T}, p::Par{T}) where T <: Real
    mean(f2.(xx, yy, p))
end

function h2(n::Int, p::Par{Float64})
    out = 0.0
    for i in 1:n
        xx = rand(500, 500)
        yy = rand(500, 500)
        out += g2(xx, yy, p)/n
    end
end


p = Par(0.1, 0.2)
