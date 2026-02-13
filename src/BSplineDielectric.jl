"""
================================================================================
Recursive B-splines to model functions that fulfill Kramers Kronig relations
Implementation of Johs and Hale, doi:10.1002/pssa.200777754, 2008
================================================================================
"""
module BSplineDielectric

export MDF, eps2basis, eps1basis, nknots, degree, nbasis, knotpos, knotgen,
    KnotBoundary, Clamped, PaddedLinear

const IdxVal = Union{Integer, AbstractVector{<:Integer}, AbstractUnitRange{<:Integer}}
const RealVal = Union{Real, AbstractVector{<:Real}}

"""
(MDF) model dielectric function
subtype of Function abstract type

### Fields
- `f`: model function
- `k`: degree of B-splines
- `t`: knots

### Methods
- `(f::MDF)(i::IdxVal, x::RealVal)`: basis vector or matrix for given knot index/indices and wavelength/wavenumber/frequency position(s) 
- `(f::MDF)(x::RealVal)`: full basis matrix (design matrix) for wavelength/wavenumber/frequency position(s)
- `knotpos(f::MDF)`: knot positions
- `degree(f::MDF)`: degree of polynomial
- `nknots(f::MDF)`: number of knots
- `nbasis(f::MDF)`: number of basis functions
where
- `IdxVal = Union{Integer, AbstractVector{<:Integer}, AbstractUnitRange{<:Integer}}`
- `RealVal = Union{Real, AbstractVector{<:Real}}`
The knots and number of knots can be retrieved with `knots` and `nknots`, respectively.
    """
struct MDF{F,K<:Integer,R<:Real} <: Function
    f::F
    k::K
    t::Vector{R}
end

(f::MDF)(i::IdxVal, x::RealVal) = f.f.(f.k, transpose(i), x)
(f::MDF)(x::RealVal) = f.f.(f.k, transpose(1:(nknots(f)-f.k-1)), x) # design matrix

knotpos(f::MDF) = f.t
degree(f::MDF) = f.k
nknots(f::MDF) = length(f.t)
nbasis(f::MDF) = nknots(f) - degree(f) - 1


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

abstract type KnotBoundary end
struct Clamped <: KnotBoundary end
struct PaddedLinear <: KnotBoundary end


"""
    function knotgen(::KnotGen, k::Integer, x::Vector{<:Real}, ...)

Function for knot generation.

### Parameters
- `::KnotBoundary`: struct which determines the method of calculation
- `k`::Integer: degree
- `x`::Vector{<:Real}: break points (interior knots)
- `...`: additional arguments

### Methods
- for `Clamped()`: augmentation results in `k+1` repeated values of the extrema of break points `x`
- for `PaddedLinear()`: pads both ends by `f * (1:k) * dx`, where `dx` at each end is estimated from the interval of the last two extrema of `x`.
- for `f::MDF`: returns knot positions for struct `MDF`

### Returns
Vector of knot positions
"""
function knotgen end

function knotgen(::Clamped, k::Integer, x::RealVal)
    a = x[begin]
    b = x[end]
    interior = x[begin+1:end-1]
    return vcat(
        fill(a, k+1), 
        interior, 
        fill(b, k+1)
    )
end

function knotgen(::PaddedLinear, k::Integer, x::RealVal, f::Real=1e2)
    a = x[begin]
    b = x[end]
    dx1 = only(diff(x[begin:begin+1]))
    dx2 = only(diff(x[end-1:end]))
    return vcat(
        a .- f .* (k:-1:1) .* dx1, 
        x, 
        b .+ f .* (1:k) .* dx2
    )
end

function isclamped(k::Int, t::AbstractVector{<:Real})
    a = t[1]
    b = t[end]
    all(t[1:k+1] .== a) && all(t[end-k:end] .== b)
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
function eps2basis(k::Integer, t::AbstractVector{<:Real})
    n = length(t)
    # B⁰(i::Integer, x::Real) = t[i] <= x < t[i+1] ? 1.0 : 0.0                  # Eq. 1
    function B⁰(i::Integer, x::Real)
        if t[i] <= x < t[i+1] return 1.0 end
        if isapprox(x, t[end]) && isapprox(x, t[i+1]) return 1.0 end
        return 0.0
    end    
    function B(k::Integer, i::Integer, x::Real)
        if (i > n - k - 1) return 0.0 end
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
    return MDF(B, k, t)
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
function eps1basis(k::Integer, t::AbstractVector{<:Real})
    n = length(t)
    function I¹(i::Integer, ω::Real)
        #= 
        Offers numerical stability in case of extremely close knot spacings,
        or `x` extremely close to or extremely far from a knot location.
        =#
        if ω > t[i+1] && t[i+1] > t[i] && t[i+2] > t[i+1]
            u₀ = (ω - t[i+1]) / (t[i+1] - t[i])                                         # Eq. 13
            u₁ = (ω - t[i+1]) / (t[i+2] - t[i+1])                                       # Eq. 13
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
        if (i > n - k - 1) return 0.0 end
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
                @rdiv((t[i+k+1] - ω) / (t[i+k+1] - t[i+1])) * I(k-1, i+1, ω)  # Eq. 11
        end
    end
    ϕ(k::Integer, i::Integer, ω::Real) = 1 / π * (I(k, i, ω) + I(k, i, -ω))   # Eq. 8
    return MDF(ϕ, k, t)
end

end # module
