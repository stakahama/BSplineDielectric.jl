"""
================================================================================
Recursive B-splines to model functions that fulfill Kramers Kronig relations
Implementation of Johs and Hale, doi:10.1002/pssa.200777754, 2008
================================================================================
"""
module BSplineDielectric

export MDF, eps2basis, eps1basis, knots, nknots, knotpos, Clamped, ExtendedEven

IdxVal = Union{Integer, AbstractVector{<:Integer}, AbstractUnitRange{<:Integer}}
RealVal = Union{Real, AbstractVector{<:Real}}

#=
Notation
| symbol            | description       |
|-------------------+-------------------|
| \(k\)             | polynomial degree |
| \(k+1\)           | order of spline   |
| \(x\), \(\omega\) | domain            |
| \(t\)             | knot              |
| \(n\)             | number of knots   |

- \(n-k-1\): number of control points; number of columns of design matrix
- \([t_{k}, t_{n-k-1}]\): domain of data
- \([t_{i}, t_{i+k+1})\): interval over which values are nonzero for each basis spline (`eps2basis`)
=#

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
- `(f::MDF)(x::RealVal)`: full basis matrix (design matrix) for wavelength/wavenumber/frequency position(s)
- `knotpos(f::MDF)`: knot positions
- `nknots(f::MDF)`: number of knots
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
(f::MDF)(x::RealVal) = f.f.(f.k, transpose(1:f.n-f.k-1), x) # design matrix

knotpos(f::MDF) = f.t
nknots(f::MDF) = f.n

"""
Conditional product operator
    ⊗(x::Real, y::Real)

`⊗` is a product operator that adheres to the rule that 0 ⋅ log(anything) = 0 (from Johs and Hale). 
"""
⊗(x::Real, y::Real) = isapprox(x, 0) ? 0.0 : x * y

"""
Ratios and logarithm of ratios from knot positions
    @rdiv((y2 - y1) / (x2 - x1))
    @rlog((y2 - y1) / (x2 - x1)) or @rlog(abs((y2 - y1) / (x2 - x1)))

### Parameters
`expr`: symbolic expression of the form `(y2 - y1) / (x2 - x1)`, or for @rlog also `abs((y2 - y1) / (x2 - x1))`

### Returns
0.0 if x1 == x2 or y2 == y1.
"""
macro rdiv(expr)
    # Extract variables
    den = expr.args[3]
    x2, x1 = den.args[2], den.args[3]

    # Build and escape the safe-if AST
    return esc(:(isapprox($(x2), $(x1)) ? 0.0 : $(expr)))
end

macro rlog(expr)
    # Check if ratio is enclosed in abs()
    if expr.args[1] === :/
        args = expr.args
    elseif expr.args[1] === :abs
        args = expr.args[2].args
    end

    # Extract variables    
    num, den = args[2], args[3]
    y2, y1 = num.args[2], num.args[3]
    x2, x1 = den.args[2], den.args[3]

    # Build and escape the safe-if AST
    return esc(:(isapprox($(y2), $(y1)) || isapprox($(x2), $(x1)) ? 0.0 : log($(expr))))
end

"""
Knot generation
    function knotpos(::KnotGen, k::Integer, x::Vector{<:Real}, ...)
    function knotpos(f::MDF)

### Parameters
- `::KnotGen`: struct which determines the method of calculation
- `k`::Integer: degree
- `x`::Vector{<:Real}: break points (interior knots)
- `...`: additional arguments

### Methods
- for `Clamped()`: augmentation results in `k+1` repeated values of the extrema of break points `x`
- for `ExtendedEven()`: extends by `f * (1:k) * dx`` on either side of `x` where `dx` is estimated from the interval of `x` (assuming evenly spaced)
- for `f::MDF`: returns knot positions for struct `MDF`

### Returns
Vector of knot positions
"""
abstract type KnotGen end
struct Clamped <: KnotGen end
struct ExtendedEven <: KnotGen end

knotpos(::Clamped, k::Integer, x::RealVal) = vcat(
    fill(x[1], k), 
    x, 
    fill(x[end], k)
)

function knotpos(::ExtendedEven, k::Integer, x::RealVal, f::Real=1e2)
    dx = only(diff(x[1:2]))
    return vcat(
        x[1] .- f .* (k:-1:1) .* dx, 
        x, 
        x[end] .+ f .* (1:k) .* dx
    )
end
                            
"""
Basis functions for ε₂ (imaginary dielectric function)

### Parameters
- `k`: degree of B-splines
- `t`: vector of knot positions

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
    n = length(t)
    B⁰(i::Integer, x::Real) = t[i] <= x < t[i+1] ? 1.0 : 0.0            # Eq. 1
    function B(k::Integer, i::Integer, x::Real)
        if (i+k+1 > n) return 0.0 end
        if k == 0
            return B⁰(i, x)
        else
            #=
            Checks added to return zero when the interval is zero,
            which happens when knots are estimated from smoothing algorithms.
            =#
            return @rdiv((x - t[i]) / (t[i+k] - t[i])) * B(k-1, i, x) +
                @rdiv((t[i+k+1] - x) / (t[i+k+1] - t[i+1])) * B(k-1, i+1, x) # Eq. 2
        end
    end
    return MDF(B, k, t, n)
end

"""
Basis functions for ε₁ (real dielectric function)

### Parameters
- `k`: degree of B-splines
- `t`: vector of knots

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
    n = length(t)
    function I¹(i::Integer, ω::Real)
        #= 
        Offers numerical stability in case of extremely close knot spacings,
        or `x` extremely close to or extremely far from a knot location.
        =#
        if ω > t[i+1] && t[i+1] > t[i] && t[i+2] > t[i+1]
            u₀ = (ω - t[i+1]) / (t[i+1] - t[i])                                        # Eq. 13
            u₁ = (ω - t[i+1]) / (t[i+2] - t[i+1])                                      # Eq. 13
            ωₛ = abs((ω - t[i+1]) / (t[i+2] - t[i]))                                    # Eq. 15
            if abs(u₀) < 1e-9 || abs(u₁) < 1e-9
                return log((t[i+2] - t[i+1]) / (t[i+1] - t[i])) -
                    (u₀ + u₁) + u₀ * log(u₀) + u₁ * log(u₁)                             # Eq. 14
            elseif ωₛ > 100
                return mapreduce(n -> ((-1 / u₀)^n - (1 / u₁)^n) / (n * (n+1)), +, 1:5) # Eq. 15
            end
        end
        #= 
        Canonical formulation.
        =#
        if isapprox(ω, t[i+1])
            return @rlog((t[i+2] - t[i+1]) / (t[i+1] - t[i]))
        end
        α₀ = @rdiv (ω - t[i]) / (t[i+1] - t[i])
        α₁ = @rdiv (t[i+2] - ω) / (t[i+2] - t[i+1])
        return α₀ ⊗ log(abs(1 - 1 / α₀)) - α₁ ⊗ log(abs(1 - 1 / α₁))
    end
    function I(k::Integer, i::Integer, ω::Real)
        if (i+k+1 > n) return 0.0 end
        if k == 1 
            return I¹(i, ω)
        else
            #=
            Checks added to return zero when the interval is zero,
            which happens when knots are estimated from smoothing algorithms.
            In contrast to the eps2 spline defined by B,
            both summation terms must be finite for results to be meaningful.
            =#
            return @rdiv((ω - t[i]) / (t[i+k] - t[i])) * I(k-1, i, ω) +
                @rdiv((t[i+k+1] - ω) / (t[i+k+1] - t[i+1])) * I(k-1, i+1, ω)      # Eq. 11
        end
    end
    ϕ(k::Integer, i::Integer, ω::Real) = 1 / π * (I(k, i, ω) + I(k, i, -ω))  # Eq. 8
    return MDF(ϕ, k, t, n)
end

end # module
