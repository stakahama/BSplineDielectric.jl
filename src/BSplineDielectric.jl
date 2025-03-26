#=
================================================================================
Recursive B-splines to model functions that fulfill Kramers Kronig relations
Implementation of Johs and Hale, doi:10.1002/pssa.200777754, 2008
================================================================================
=#

module BSplineDielectric

export MDF, expand, eps2basis, eps1basis

#=
(MDF) model dielectric function
also includes information on number of knots
=#
struct MDF
    basis::Function
    n::Integer
end

function expand(x::MDF)
    return x.basis(1:x.n)
end

"""
Basis functions for ε₂ (imaginary dielectric function)

### Parameters
k: degree of B-splines
t: vector of knots
x: wavenumber or frequency

### Return value
struct of type MDF

### Example
B = eps2basis(k, i, x)

## at first knot
plot(x, B.basis(1))

## at first two knots
plot(x, B.basis(1:2))

## full basis expansion
Bmatrix = expand(B)
plot(x, Bmatrix * ones(size(Bmatrix, 2)))
"""
function eps2basis(k, t, x)
    Nₜ = length(t)
    B⁰(i, x) = t[i] <= x < t[i+1] ? 1 : 0
    function B(k, i, x)
        if (i+k+1 > Nₜ); return 0; end
        if k == 0
            B⁰(i, x)
        else
            (x - t[i]) / (t[i+k] - t[i]) * B(k-1, i, x) +
                (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * B(k-1, i+1, x)
        end
    end
    # return B.(k, transpose(1:Nₜ), x)
    return MDF(i -> B.(k, transpose(i), x), Nₜ)
end

"""
Basis functions for ε₁ (real dielectric function)

### Parameters
k: degree of B-splines
t: vector of knots
x: wavenumbers

### Return value
struct of type MDF

### Example
ϕ = eps1basis(k, i, x)

## at first knot
plot(x, ϕ.basis(1))

## at first two knots
plot(x, ϕ.basis(1:2))

## full basis expansion
ϕmatrix = expand(ϕ)
plot(x, ϕmatrix * ones(size(ϕmatrix, 2)))
"""
function eps1basis(k, t, x)
    Nₜ = length(t)
    xlogx(x) = isapprox(x, 0) ? 0 : x * log(x)
    xlogy(x, y) = isapprox(x, 0) ? 0 : x * log(y)
    function I¹(i, ω)
        u₀ = (ω - t[i+1]) / (t[i+1] - t[i])
        u₁ = (ω - t[i+1]) / (t[i+2] - t[i+1])
        ωₛ = (ω - t[i+1]) / (t[i+2] - t[i])
        if abs(u₀) < 1e-9 || abs(u₁) < 1e-9
            log((t[i+2] - t[i+1]) / (t[i+1] - t[i])) -
                (u₀ + u₁) + xlogx(u₀) + xlogx(u₁)            
        elseif abs(ωₛ) > 100
            mapreduce(n -> ((-1 / u₀)^n - (1 / u₁)^n) / n(n+1), +, 1:5)
        else
            if isapprox(ω, t[i+1])
                return log((t[i+2] - t[i+1]) / (t[i+1] - t[i]))
            end
            α₀ = (ω - t[i]) / (t[i+1] - t[i])
            α₁ = (t[i+2] - ω) / (t[i+2] - t[i+1])
            xlogy(α₀, abs(1 - 1 / α₀)) - xlogy(α₁, abs(1 - 1 / α₁))            
        end
    end
    function I(k, i, ω)
        if (i+k+1 > Nₜ); return 0; end
        if k == 1 
            I¹(i, ω)
        else
            (ω - t[i]) / (t[i+k] - t[i]) * I(k-1, i, ω) +
                (t[i+k+1] - ω) / (t[i+k+1] - t[i+1]) * I(k-1, i+1, ω)
        end
    end
    ϕ(k, i, ω) = 1 / π * (I(k, i, ω) + I(k, i, -ω))
    # return ϕ.(k, transpose(1:Nₜ), x)
    return MDF(i -> ϕ.(k, transpose(i), x), Nₜ)
end

end
