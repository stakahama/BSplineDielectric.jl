"""
================================================================================
Recursive B-splines to model functions that fulfill Kramers Kronig relations
Implementation of Johs and Hale, doi:10.1002/pssa.200777754, 2008
================================================================================
"""
module BSplineDielectric

export MDF, eps2basis, eps1basis, knots, nknots

IdxType = Union{Integer, AbstractVector{<:Integer}, AbstractUnitRange{<:Integer}}
VarType = Union{Real, AbstractVector{<:Real}}

"""
(MDF) model dielectric function
subtype of Function abstract type

### Fields
- `f`: model function
- `k`: degree of B-splines
- `t`: knots
- `n`: number of knots

### Methods
- `(f::MDF)(i::IdxType, x::VarType)`: basis vector or matrix for given knot index/indices and wavelength/wavenumber/frequency position(s) 
- `(f::MDF)(x::VarType)`: full basis matrix for wavelength/wavenumber/frequency position(s)
where
- `IdxType = Union{Integer, AbstractVector{<:Integer}, AbstractUnitRange{<:Integer}}`
- `VarType = Union{Real, AbstractVector{<:Real}}`
The knots and number of knots can be retrieved with `knots` and `nknots`, respectively.
"""
struct MDF <: Function
    f::Function
    k::Integer
    t::Vector{<:Real}
    n::Integer
end

(f::MDF)(i::IdxType, x::VarType) = f.f.(f.k, transpose(i), x)
(f::MDF)(x::VarType) = f.f.(f.k, transpose(1:f.n), x)
knots(f::MDF) = f.t
nknots(f::MDF) = f.n

"""
Basis functions for ε₂ (imaginary dielectric function)

### Parameters
- k: degree of B-splines
- t: vector of knot positions

### Returns
struct of type MDF - can be called as a function.

### Example
    B = eps2basis(k, knots)

    ## at first knot
    plot(x, B(1, x))

    ## at first two knots
    plot(x, B(1:2, x))

    ## full basis expansion
    plot(x, B(x) * ones(nknots(B)))
"""
function eps2basis(k, t)
    N = length(t)
    B⁰(i, x) = t[i] <= x < t[i+1] ? 1 : 0
    function B(k, i, x)
        if (i+k+1 > N) return 0 end
        if k == 0
            B⁰(i, x)
        else
            (x - t[i]) / (t[i+k] - t[i]) * B(k-1, i, x) +
                (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * B(k-1, i+1, x)
        end
    end
    return MDF(B, k, t, N)
end

"""
Basis functions for ε₁ (real dielectric function)

### Parameters
- k: degree of B-splines
- t: vector of knots

### Returns
struct of type MDF - can be called as a function.

### Example
    ϕ = eps1basis(k, knots)

    ## at first knot
    plot(x, ϕ(1, x))

    ## at first two knots
    plot(x, ϕ(1:2, x))

    ## full basis expansion
    plot(x, ϕ(x) * ones(nknots(ϕ)))
"""
function eps1basis(k, t)
    N = length(t)
    xlogy(x, y = x) = isapprox(x, 0) ? 0 : x * log(y)
    function I¹(i, ω)
        u₀ = (ω - t[i+1]) / (t[i+1] - t[i])
        u₁ = (ω - t[i+1]) / (t[i+2] - t[i+1])
        ωₛ = (ω - t[i+1]) / (t[i+2] - t[i])
        if abs(u₀) < 1e-9 || abs(u₁) < 1e-9
            log((t[i+2] - t[i+1]) / (t[i+1] - t[i])) -
                (u₀ + u₁) + xlogy(u₀) + xlogy(u₁)            
        elseif abs(ωₛ) > 100
            mapreduce(n -> ((-1 / u₀)^n - (1 / u₁)^n) / (n * (n+1)), +, 1:5)
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
        if (i+k+1 > N) return 0 end
        if k == 1 
            I¹(i, ω)
        else
            (ω - t[i]) / (t[i+k] - t[i]) * I(k-1, i, ω) +
                (t[i+k+1] - ω) / (t[i+k+1] - t[i+1]) * I(k-1, i+1, ω)
        end
    end
    ϕ(k, i, ω) = 1 / π * (I(k, i, ω) + I(k, i, -ω))
    return MDF(ϕ, k, t, N)
end

end
