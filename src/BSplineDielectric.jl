"""
================================================================================
Recursive B-splines to model functions that fulfill Kramers Kronig relations
Implementation of Johs and Hale, doi:10.1002/pssa.200777754, 2008
================================================================================
"""
module BSplineDielectric

export MDF, eps2basis, eps1basis, knots, nknots

IdxVal = Union{Integer, AbstractVector{<:Integer}, AbstractUnitRange{<:Integer}}
RealVal = Union{Real, AbstractVector{<:Real}}

"""
(MDF) model dielectric function
subtype of Function abstract type

### Fields
- `f`: model function
- `k`: degree of B-splines
- `t`: knots
- `n`: number of knots

### Methods
- `(f::MDF)(i::IdxVal, x::RealVal)`: basis vector or matrix for given knot index/indices and wavelength/wavenumber/frequency position(s) 
- `(f::MDF)(x::RealVal)`: full basis matrix for wavelength/wavenumber/frequency position(s)
where
- `IdxVal = Union{Integer, AbstractVector{<:Integer}, AbstractUnitRange{<:Integer}}`
- `RealVal = Union{Real, AbstractVector{<:Real}}`
The knots and number of knots can be retrieved with `knots` and `nknots`, respectively.
"""
struct MDF <: Function
    f::Function
    k::Integer
    t::Vector{<:Real}
    n::Integer
end

(f::MDF)(i::IdxVal, x::RealVal) = f.f.(f.k, transpose(i), x)
(f::MDF)(x::RealVal) = f.f.(f.k, transpose(1:f.n), x)
knots(f::MDF) = f.t
nknots(f::MDF) = f.n

# ⊘(x::Real, y::Real) = isapprox(y, 0) ? 0 : x / y
⊗(x::Real, y::Real) = isapprox(x, 0) ? 0 : x * y
⊕(x::Real, y::Real) = isfinite(x) && isfinite(y) ? x + y : 0
# ⊛(x::Real, y::Real) = isfinite(x) ? x * y : 0
    
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
function eps2basis(k::Integer, t::RealVal)
    N = length(t)
    B⁰(i::Integer, x::Real) = t[i] <= x < t[i+1] ? 1 : 0                # Eq. 1
    function B(k::Integer, i::Integer, x::Real)
        if (i+k+1 > N) return 0 end
        if k == 0
            B⁰(i, x)
        else
            #=
            Checks added to return zero when the interval is zero,
            which happens when knots are estimated from smoothing algorithms.
            =#
            (isapprox(t[i+k], t[i]) ? 0 : (x - t[i]) / (t[i+k] - t[i]) * B(k-1, i, x)) +
               (isapprox(t[i+k+1], t[i+1]) ? 0 : (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * B(k-1, i+1, x)) # Eq. 2
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
function eps1basis(k::Integer, t::RealVal)
    N = length(t)
    function I¹(i::Integer, ω::Real)
        u₀ = (ω - t[i+1]) / (t[i+1] - t[i])                                  # Eq. 13
        u₁ = (ω - t[i+1]) / (t[i+2] - t[i+1])                                # Eq. 13
        ωₛ = abs((ω - t[i+1]) / (t[i+2] - t[i]))                             # Eq. 15
        if abs(u₀) < 1e-9 || abs(u₁) < 1e-9
            log((t[i+2] - t[i+1]) / (t[i+1] - t[i])) -
                (u₀ + u₁) + u₀ ⊗ log(u₀) + u₁ ⊗ log(u₁)                      # Eq. 14
        elseif ωₛ > 100
            mapreduce(n -> ((-1 / u₀)^n - (1 / u₁)^n) / (n * (n+1)), +, 1:5) # Eq. 15
        else                                                                 # Eq. 12
            if isapprox(ω, t[i+1])
                return log((t[i+2] - t[i+1]) / (t[i+1] - t[i]))
            end
            α₀ = (ω - t[i]) / (t[i+1] - t[i])
            α₁ = (t[i+2] - ω) / (t[i+2] - t[i+1])
            α₀ ⊗ log(abs(1 - 1 / α₀)) - α₁ ⊗ log(abs(1 - 1 / α₁))
        end
    end
    function I(k::Integer, i::Integer, ω::Real)
        if (i+k+1 > N) return 0 end
        if k == 1 
            I¹(i, ω)
        else
            #=
            Checks added to return zero when the interval is zero,
            which happens when knots are estimated from smoothing algorithms.
            In contrast to the eps2 spline defined by B,
            both terms must be finite for results to be meaningful.
            This is not the most efficient implementation as it checks only
            after the fact.
            =#
            (ω - t[i]) / (t[i+k] - t[i]) * I(k-1, i, ω) ⊕
                (t[i+k+1] - ω) / (t[i+k+1] - t[i+1]) * I(k-1, i+1, ω)        # Eq. 11
        end
    end
    ϕ(k::Integer, i::Integer, ω::Real) = 1 / π * (I(k, i, ω) + I(k, i, -ω))  # Eq. 8
    return MDF(ϕ, k, t, N)
end

end
