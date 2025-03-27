BSplineDielectric.jl
==============

This package implements recursive B-splines that fulfill Kramers Kronig relations as described by Johs and Hale (doi:[10.1002/pssa.200777754](https://doi.org/10.1002/pssa.200777754), 2008). 

Johs and Hale illustrated their use for modeling dielectric functions, but they can be used to represent any function that fulfills the Kramers-Kronig relation - such as the complex refractive index and Clausius-Mossotti factor.

Installation.
```julia
using Pkg
Pkg.add(url = "https://github.com/stakahama/BSplineDielectric.jl")
```

Example usage.
```julia
## Libraries
using Plots
using BSplineDielectric

## Define x values and knot positions; generate curve
function triangle(x)
    2 < x ≤ 3 ? x - 2 :
    3 < x ≤ 4 ? 4 - x :
    0
end

x = collect(range(0, 6, length = 101))
iknots = range(1, length(x), step = 2)
curve = triangle.(x)

## Identify knot positions and generate spline bases
B = eps2basis(1, x[iknots])
ϕ = eps1basis(1, x[iknots])

## View first terms of spline basis
plot(x, B(1, x), legend = false)
plot(x, B(1:2, x), legend = false)

## Get number of knots (maximum allowable index)
nknots(B)  # same with nknots(ϕ)

## Get knot positions
knots(B)  # same with knots(ϕ)

## Get full basis matrices evaluated at x
Bmat = B(x)
ϕmat = ϕ(x)

## Estimate spline coefficients
coef =  Bmat \ curve

## Plot fits
plot(x, curve, label = "ε₂")
plot!(x, Bmat * coef, label = "ε₂ model")
plot!(x, ϕmat * coef, label = "ε₁ model")

## Apply same coefficients on B-splines of degree 3
B3 = eps2basis(3, x[iknots])
ϕ3 = eps1basis(3, x[iknots])

## Plot
plot(x, B3(x) * coef, label = "ε₂ model")
plot!(x, ϕ3(x) * coef, label = "ε₁ model")
```
