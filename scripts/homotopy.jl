## homotopy.jl --- 

###_____________________________________________________________________________
##
##  Load Packages
###_____________________________________________________________________________

import Pkg
Pkg.add("HomotopyContinuation")
Pkg.add("DelimitedFiles")
using HomotopyContinuation, LinearAlgebra, DelimitedFiles


# This script was tested with Julia v1.5.3.

###_____________________________________________________________________________
## 
##  Define variables and some constraints
###_____________________________________________________________________________

@var  θ[1:5] 
@var f₋₋₋₋ f₋₋₋₊ f₋₋₊₋ f₋₋₊₊ f₋₊₋₋ f₋₊₋₊ f₋₊₊₋ f₋₊₊₊ f₊₋₋₋ f₊₋₋₊ f₊₋₊₋ f₊₋₊₊ f₊₊₋₋ f₊₊₋₊ f₊₊₊₋ f₊₊₊₊  θ₁ θ₂ θ₃ θ₄ θ₅ 

f₋₋₋₋ = (1/16) * (1 + θ₁*θ₅*θ₂ + θ₁*θ₃ + θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ + θ₂*θ₄ + θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)
f₋₋₋₊ = (1/16) * (1 + θ₁*θ₅*θ₂ + θ₁*θ₃ - θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ - θ₂*θ₄ - θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₋₋₊₋ = (1/16) * (1 + θ₁*θ₅*θ₂ - θ₁*θ₃ + θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ + θ₂*θ₄ - θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₋₋₊₊ = (1/16) * (1 + θ₁*θ₅*θ₂ - θ₁*θ₃ - θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ - θ₂*θ₄ + θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)
f₋₊₋₋ = (1/16) * (1 - θ₁*θ₅*θ₂ + θ₁*θ₃ + θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ - θ₂*θ₄ + θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₋₊₋₊ = (1/16) * (1 - θ₁*θ₅*θ₂ + θ₁*θ₃ - θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ + θ₂*θ₄ - θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)
f₋₊₊₋ = (1/16) * (1 - θ₁*θ₅*θ₂ - θ₁*θ₃ + θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ - θ₂*θ₄ - θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)
f₋₊₊₊ = (1/16) * (1 - θ₁*θ₅*θ₂ - θ₁*θ₃ - θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ + θ₂*θ₄ + θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₊₋₋₋ = (1/16) * (1 - θ₁*θ₅*θ₂ - θ₁*θ₃ - θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ + θ₂*θ₄ + θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₊₋₋₊ = (1/16) * (1 - θ₁*θ₅*θ₂ - θ₁*θ₃ + θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ - θ₂*θ₄ - θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)
f₊₋₊₋ = (1/16) * (1 - θ₁*θ₅*θ₂ + θ₁*θ₃ - θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ + θ₂*θ₄ - θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)
f₊₋₊₊ = (1/16) * (1 - θ₁*θ₅*θ₂ + θ₁*θ₃ + θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ - θ₂*θ₄ + θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₊₊₋₋ = (1/16) * (1 + θ₁*θ₅*θ₂ - θ₁*θ₃ - θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ - θ₂*θ₄ + θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)
f₊₊₋₊ = (1/16) * (1 + θ₁*θ₅*θ₂ - θ₁*θ₃ + θ₁*θ₅*θ₄ - θ₂*θ₅*θ₃ + θ₂*θ₄ - θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₊₊₊₋ = (1/16) * (1 + θ₁*θ₅*θ₂ + θ₁*θ₃ - θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ - θ₂*θ₄ - θ₃*θ₅*θ₄ - θ₁*θ₃*θ₂*θ₄)
f₊₊₊₊ = (1/16) * (1 + θ₁*θ₅*θ₂ + θ₁*θ₃ + θ₁*θ₅*θ₄ + θ₂*θ₅*θ₃ + θ₂*θ₄ + θ₃*θ₅*θ₄ + θ₁*θ₃*θ₂*θ₄)



f₋₋₋₋ = (1/16) * (1 + τ)
f₋₋₋₊ = (1/16) * (1 + τ)
f₋₋₊₋ = (1/16) * (1 + τ)
f₋₋₊₊ = (1/16) * (1 + τ)
f₊₊₋₋ = (1/16) * (1 + τ)
f₊₊₋₊ = (1/16) * (1 + τ)
f₊₊₊₋ = (1/16) * (1 + τ)
f₊₊₊₊ = (1/16) * (1 + τ)


f₋₊₋₋ = (1/16) * (1 - τ)
f₋₊₋₊ = (1/16) * (1 - τ)
f₋₊₊₋ = (1/16) * (1 - τ)
f₋₊₊₊ = (1/16) * (1 - τ)
f₊₋₋₋ = (1/16) * (1 - τ)
f₊₋₋₊ = (1/16) * (1 - τ)
f₊₋₊₋ = (1/16) * (1 - τ)
f₊₋₊₊ = (1/16) * (1 - τ)


N₊= sum([n[i] for i in [1 2 3 4 13 14 15 16]])
N₋= sum([n[i] for i in 5:12])
@var τ λ
log_L = log((1+τ)/16)*N₊+log((1-τ)/16)*N₋
Lagrangian=log_L + λ*


f=[f₋₋₋₋, f₋₋₋₊, f₋₋₊₋, f₋₋₊₊, f₋₊₋₋, f₋₊₋₊, f₋₊₊₋, f₋₊₊₊, f₊₋₋₋, f₊₋₋₊, f₊₋₊₋,f₊₋₊₊, f₊₊₋₋, f₊₊₋₊, f₊₊₊₋, f₊₊₊₊]
@var f[1:16]
f[1]
f[2]

## S = sum(f)-1
T = θ₁*θ₅*θ₂ - 1/2
Z = [θ₃; θ₄]
H = f-h


###_____________________________________________________________________________
##
## Define the critical equations + constraints
###_____________________________________________________________________________

@var n[1:16]
@var logL Lagrangian λ[1:4]
logL = sum(n.*log.(f))
logL
Lagrangian = logL + λ'*[S;T;Z]

dL=differentiate(Lagrangian, θ)

dL
F = Diagonal(f)

corrected_dL = F*dL
# system of equations

C = System([corrected_dL; S; T; H; Z], variables = [f; θ; λ], parameters = n)
n₀ = rand(Float64, 16)

res = solve(C; target_parameters = n₀) 

###_____________________________________________________________________________
##
##  Some new constraints
###_____________________________________________________________________________

# These constraints come from Example 2.6 in Ben's Paper. The following two
# constraints apply in the case that the tree has topology 12|34
(f[1] - f[2] - f[3] + f[4] - f[5] + f[6] + f[7] - f[8] - f[9] + f[10] + f[11] - f[12] + f[13] - f[14] - f[15] + f[16])*(f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16]) - (f[1] - f[2] - f[3] + f[4] + f[5] - f[6] - f[7] + f[8] + f[9] - f[10] - f[11] + f[12] + f[13] - f[14] - f[15] + f[16])*(f[1] + f[2] + f[3] + f[4] - f[5] - f[6] - f[7] - f[8] - f[9] - f[10] - f[11] - f[12] + f[13] + f[14] + f[15] + f[16])


(f[1] + f[2] - f[3] - f[4] + f[5] + f[6] - f[7] - f[8] - f[9] - f[10] + f[11] + f[12] - f[13] - f[14] + f[15] + f[16])*(f[1] - f[2] + f[3] - f[4] - f[5] + f[6] - f[7] + f[8] + f[9] - f[10] + f[11] - f[12] - f[13] + f[14] - f[15] + f[16]) - (f[1] + f[2] - f[3] - f[4] - f[5] - f[6] + f[7] + f[8] + f[9] + f[10] - f[11] - f[12] - f[13] - f[14] + f[15] + f[16])*(f[1] - f[2] + f[3] - f[4] + f[5] - f[6] + f[7] - f[8] - f[9] + f[10] - f[11] + f[12] - f[13] + f[14] - f[15] + f[16])



(f₁ - f₂ - f₃ + f₄ - f₅ + f₆ + f₇ - f₈ - f₉ + f₁₀ + f₁₁ - f₁₂ + f₁₃ - f₁₄ - f₁₅ + f₁₆)*(f₁ + f₂ + f₃ + f₄ + f₅ + f₆ + f₇ + f₈ + f₉ + f₁₀ + f₁₁ + f₁₂ + f₁₃ + f₁₄ + f₁₅ + f₁₆) - (f₁ - f₂ - f₃ + f₄ + f₅ - f₆ - f₇ + f₈ + f₉ - f₁₀ - f₁₁ + f₁₂ + f₁₃ - f₁₄ - f₁₅ + f₁₆)*(f₁ + f₂ + f₃ + f₄ - f₅ - f₆ - f₇ - f₈ - f₉ - f₁₀ - f₁₁ - f₁₂ + f₁₃ + f₁₄ + f₁₅ + f₁₆)


(f₁ + f₂ - f₃ - f₄ + f₅ + f₆ - f₇ - f₈ - f₉ - f₁₀ + f₁₁ + f₁₂ - f₁₃ - f₁₄ + f₁₅ + f₁₆)*(f₁ - f₂ + f₃ - f₄ - f₅ + f₆ - f₇ + f₈ + f₉ - f₁₀ + f₁₁ - f₁₂ - f₁₃ + f₁₄ - f₁₅ + f₁₆) - (f₁ + f₂ - f₃ - f₄ - f₅ - f₆ + f₇ + f₈ + f₉ + f₁₀ - f₁₁ - f₁₂ - f₁₃ - f₁₄ + f₁₅ + f₁₆)*(f₁ - f₂ + f₃ - f₄ + f₅ - f₆ + f₇ - f₈ - f₉ + f₁₀ - f₁₁ + f₁₂ - f₁₃ + f₁₄ - f₁₅ + f₁₆)


# The next two constraints are for the case when the true topology is 13|24
