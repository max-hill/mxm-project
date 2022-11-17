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

# @var τ

# f₋₋₋₋ = (1/16) * (1 + τ)
# f₋₋₋₊ = (1/16) * (1 + τ)
# f₋₋₊₋ = (1/16) * (1 + τ)
# f₋₋₊₊ = (1/16) * (1 + τ)
# f₊₊₋₋ = (1/16) * (1 + τ)
# f₊₊₋₊ = (1/16) * (1 + τ)
# f₊₊₊₋ = (1/16) * (1 + τ)
# f₊₊₊₊ = (1/16) * (1 + τ)


# f₋₊₋₋ = (1/16) * (1 - τ)
# f₋₊₋₊ = (1/16) * (1 - τ)
# f₋₊₊₋ = (1/16) * (1 - τ)
# f₋₊₊₊ = (1/16) * (1 - τ)
# f₊₋₋₋ = (1/16) * (1 - τ)
# f₊₋₋₊ = (1/16) * (1 - τ)
# f₊₋₊₋ = (1/16) * (1 - τ)
# f₊₋₊₊ = (1/16) * (1 - τ)


# N₊= sum([n[i] for i in [1 2 3 4 13 14 15 16]])
# N₋= sum([n[i] for i in 5:12])
# @var τ λ
# log_L = log((1+τ)/16)*N₊+log((1-τ)/16)*N₋
# Lagrangian=log_L + λ*


# define the vector of probabilities
f=[f₋₋₋₋, f₋₋₋₊, f₋₋₊₋, f₋₋₊₊, f₋₊₋₋, f₋₊₋₊, f₋₊₊₋, f₋₊₊₊, f₊₋₋₋, f₊₋₋₊, f₊₋₊₋,f₊₋₊₊, f₊₊₋₋, f₊₊₋₊, f₊₊₊₋, f₊₊₊₊]
@var f[1:16]


# constraints 
g = sum(f) - 1
C1 = (f[1] - f[2] - f[3] + f[4] - f[5] + f[6] + f[7] - f[8] - f[9] + f[10] + f[11] - f[12] + f[13] - f[14] - f[15] + f[16])*(f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16]) - (f[1] - f[2] - f[3] + f[4] + f[5] - f[6] - f[7] + f[8] + f[9] - f[10] - f[11] + f[12] + f[13] - f[14] - f[15] + f[16])*(f[1] + f[2] + f[3] + f[4] - f[5] - f[6] - f[7] - f[8] - f[9] - f[10] - f[11] - f[12] + f[13] + f[14] + f[15] + f[16])

C2 = (f[1] + f[2] - f[3] - f[4] + f[5] + f[6] - f[7] - f[8] - f[9] - f[10] + f[11] + f[12] - f[13] - f[14] + f[15] + f[16])*(f[1] - f[2] + f[3] - f[4] - f[5] + f[6] - f[7] + f[8] + f[9] - f[10] + f[11] - f[12] - f[13] + f[14] - f[15] + f[16]) - (f[1] + f[2] - f[3] - f[4] - f[5] - f[6] + f[7] + f[8] + f[9] + f[10] - f[11] - f[12] - f[13] - f[14] + f[15] + f[16])*(f[1] - f[2] + f[3] - f[4] + f[5] - f[6] + f[7] - f[8] - f[9] + f[10] - f[11] + f[12] - f[13] + f[14] - f[15] + f[16])

C = [g;C1;C2] # a vector of all the constraints

# likelihood 
@var n[1:16]
@var logL Lagrangian λ[1:3]
logL = sum(n.*log.(f))
logL
Lagrangian = logL + λ'*[g;C1;C2]

dL = dL=differentiate(Lagrangian, f)

F = Diagonal(f)

corrected_dL = F*dL


# system of equations -- takes 1.5 hours to solve

C = System([corrected_dL; g; C1; C2], variables = [f; λ], parameters = n)
n₀ = rand(Float64, 16)

res = solve(C; target_parameters = n₀) 

real_solutions(res)



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


# tbd


###_____________________________________________________________________________
##
##  Critical Point Classification
###_____________________________________________________________________________

J = differentiate(C,f) # where C is a vector of the constraints (e.g.
                       # C=[constraint1; constraint2; constraint3]). J is the
                       # jacobian of the constraints
m = 3 # number of constraints (3 in our case)
N = 4 # number of leaves (4 in our case)
H = differentiate(differentiate(Lagrangian,f),f) # The Hessian of the Lagrangian.

random_n = n₀ # the data
n

###
function classify_critical_point(critical_point::Array)

    "This function implements the bordered determinental criterion (i.e. 2nd
    derivative test) to test if a given real critical point is a local maximum.
    Specifically, we follow Theorem 12 on page 155 of Magnus and Neudecker, Jan
    R. and Heinz. Matrix Differential Calculus with Applications in Statistics
    and Econometrics. Second (revised) edition, 1999. It returns 1 if the point
    is a maximum, -1 if the point is a minimum, and 0 if the test is
    indeterminate."

    c=critical_point[1:2^N] # The pᵢ's
    d=critical_point[(2^N+1):2^N+m] # The λᵢ's
    M=fill(0,m,m) # Square m×m matrix of zeros.

    # First, we need to test that the the m×m upper-left square submatrix of J
    # has nonzero determinant when evaluated at the critical point. Magnus et al
    # claim this condition is necessary but that we can always relabel our
    # variables so that this condition holds. I have not implemented relabeling
    # but will if it becomes necessary to do so. For now, the function just
    # tests to see if the condition holds and gives up if it does not.

    Jₘ=view(J,:,1:m)
    Jₘ=subs(Jₘ,[f;λ;n]=>[c;d;random_n]) # evaluate it at the critical point
    Jₘ=convert(Array{Float32,2},Jₘ) # convert it to float
    if det(Jₘ)==0
        return "Error: unmet condition. Indeterminate."
    else
        tmp_vector=[]
        for r in (m+1):2^N
            Hᵣ=view(H,1:r,1:r) # r×r submatrix of the Hessian of the Lagrangian
            Jᵣ=view(J,:,1:r)   # First r columns of the Jacobian of constraints

            # Next, we construct Δ, which is a symmetric (m+r)×(m+r) minor of
            # the (m+N)×(m+N) bordered Hessian. (In particular, Δ is the matrix
            # defined by equation (6) in Magnus et al. on page 155.) We then
            # evaluate Δ at the critical point, and with the randomly-chosen
            # parameters used when constructing our system. For computational
            # reasons we then tell Julia to treat the entries of Δ as numbers,
            # rather than symbolic expressions, and compute its determinant.
            
            Δ_matrix = [M Jᵣ;Jᵣ' Hᵣ]
            Δ_matrix = subs(Δ_matrix, [f;λ;n]=>[c;d;random_n])
            Δ_matrix = convert(Array{Float32,2},Δ_matrix)
            Δ_det = prod(eigvals(Δ_matrix))
            
            # Collect the determinant into a temporary vector.
            tmp_vector=[tmp_vector;Δ_det]
        end
        
        # Finally, check the determinental condition for classification.
        if [tmp_vector[i]*(-1)^(i+m)>0 for i in 1:(2^N-m)] == fill(true, 2^N-m)
            return 1 # yes, maximum. note i=r-m
        elseif [tmp_vector[i]*(-1)^m>0 for i in 1:(2^N-m)] == fill(true, 2^N-m)
            return -1 # no, minimum
        else
            return 0 # saddle-point or degenerate extrema
        end
    end
end
###
classify_critical_point(real_sols[1])
classify_critical_point(real_sols[2])
subs(logL, [n;f] => [n₀; real_sols[2][1:16]])



###_____________________________________________________________________________
## Homotopy Method for Least Square Code - Case 1
###_____________________________________________________________________________

# minimize equation on slide 12
square(n) = n * n
@var θ[1:5]
@var z[1:16]
@var α[3:4]
@var f₋₋₋₋ f₋₋₋₊ f₋₋₊₋ f₋₋₊₊ f₋₊₋₋ f₋₊₋₊ f₋₊₊₋ f₋₊₊₊ f₊₋₋₋ f₊₋₋₊ f₊₋₊₋ f₊₋₊₊ f₊₊₋₋ f₊₊₋₊ f₊₊₊₋ f₊₊₊₊  θ₁ θ₂ θ₅ α₃ α₄
θ = [θ₁, θ₂, θ₅, α₃, α₄]
τ = (1/2)

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

f=[f₋₋₋₋, f₋₋₋₊, f₋₋₊₋, f₋₋₊₊, f₋₊₋₋, f₋₊₋₊, f₋₊₊₋, f₋₊₊₊, f₊₋₋₋, f₊₋₋₊, f₊₋₊₋,f₊₋₊₊, f₊₊₋₋, f₊₊₋₊, f₊₊₊₋, f₊₊₊₊]


@var ϵ₋₋₋₋ ϵ₋₋₋₊ ϵ₋₋₊₋ ϵ₋₋₊₊ ϵ₋₊₋₋ ϵ₋₊₋₊ ϵ₋₊₊₋ ϵ₋₊₊₊ ϵ₊₋₋₋ ϵ₊₋₋₊ ϵ₊₋₊₋ ϵ₊₋₊₊ ϵ₊₊₋₋ ϵ₊₊₋₊ ϵ₊₊₊₋ ϵ₊₊₊₊ α₃ α₄

ϵ₋₋₋₋ =   θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₋₊ =   θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₊₋ = - θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₋₊₊ = - θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₋ =   θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₊ =   θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₊₋ = - θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₊₊₊ = - θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₋ = - θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₊ = - θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₊₋ =   θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₋₊₊ =   θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₋ = - θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₊ = - θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₊₋ =   θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₊₊₊ =   θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]

#define constraints: τ = 1/2
g = θ₁ *θ₂ *θ₅ - 0.5
h = [θ₁- .7633788129234264;  α₃]  # boundary cases


# Homotopy Implementation

@var Lagrangian λ[1:3] L
f_inv = map(inv,f)
L = (1/2)*sum((f_inv).*square.(z- ϵ))
#Lagrangian method 
Lagrangian = L + λ'*[g; h]
#Jacobian Matrix 
J = differentiate(Lagrangian,θ)
C = System([J; g; h], variables = [θ; λ], parameters = z)


# @z_0: input data for z 

z₀ = [-0.08984943, -0.23977824, 0.34743529, -0.33314363, -0.01655752, -0.07218132,
-0.17829711, -0.16338844, -0.3101673, -0.25123448, 0.2592691, -0.09633576,
 0.47430885, 0.38384785, -0.15680448, 0.44287663]


res = solve(C; target_parameters = z₀)

sols=real_solutions(res)
res[8]



# finding max/min value
r = real_solutions(res)

θ₁, θ₂, θ₅, α₃, α₄, λ = r[1]  # change based on which solution gives a sensible result

ϵ₋₋₋₋ = θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₋₊ = θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₊₋ = - θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₋₊₊ = - θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₋ = θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₊ = θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₊₋ = - θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₊₊₊ = - θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₋ = - θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₊ = - θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₊₋ = θ₁*α₃ + θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₋₊₊ = θ₁*α₃ - θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₋ = - θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₊ = - θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₊₋ = θ₁*α₃ - θ₂*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₊₊₊ = θ₁*α₃ + θ₂*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]


(1/2)*sum((f_inv).*square.(z₀ - ϵ))


###_____________________________________________________________________________
## Homotopy Method for Least Square Code - Case 3
###_____________________________________________________________________________

# minimize equation on slide 12
square(n) = n * n
@var θ[1:5]
@var z[1:16]
@var α[3:4]
@var f₋₋₋₋ f₋₋₋₊ f₋₋₊₋ f₋₋₊₊ f₋₊₋₋ f₋₊₋₊ f₋₊₊₋ f₋₊₊₊ f₊₋₋₋ f₊₋₋₊ f₊₋₊₋ f₊₋₊₊ f₊₊₋₋ f₊₊₋₊ f₊₊₊₋ f₊₊₊₊  θ₁ θ₂ θ₅ α₃ α₄
θ = [θ₁, θ₂, θ₅, α₃, α₄]
τ = (1/2)

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

f=[f₋₋₋₋, f₋₋₋₊, f₋₋₊₋, f₋₋₊₊, f₋₊₋₋, f₋₊₋₊, f₋₊₊₋, f₋₊₊₊, f₊₋₋₋, f₊₋₋₊, f₊₋₊₋,f₊₋₊₊, f₊₊₋₋, f₊₊₋₊, f₊₊₊₋, f₊₊₊₊]


@var ϵ₋₋₋₋ ϵ₋₋₋₊ ϵ₋₋₊₋ ϵ₋₋₊₊ ϵ₋₊₋₋ ϵ₋₊₋₊ ϵ₋₊₊₋ ϵ₋₊₊₊ ϵ₊₋₋₋ ϵ₊₋₋₊ ϵ₊₋₊₋ ϵ₊₋₊₊ ϵ₊₊₋₋ ϵ₊₊₋₊ ϵ₊₊₊₋ ϵ₊₊₊₊ α₃ α₄

ϵ₋₋₋₋ = θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₋₋₊ = - θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₋₋₊₋ = θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₋₊₊ = - θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₋₊₋₋ = θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₊₋₊ = - θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₋₊₊₋ = θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₊₊₊ = - θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₋₋₋ = - θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₋₋₊ = θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₊₋₊₋ = - θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₋₊₊ = θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₊₊₋₋ = - θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₊₋₊ = θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₊₊₊₋ = - θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₊₊₊ = θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]

#define constraints: τ = 1/2
g = θ₁ *θ₂ *θ₅ - 0.5



# Homotopy Implementation

@var Lagrangian λ L
f_inv = map(inv,f)
L = (1/2)*sum((f_inv).*square.(z- ϵ))
#Lagrangian method 
Lagrangian = L + λ*g
#Jacobian Matrix 
J = differentiate(Lagrangian,θ)
C = System([J; g], variables = [θ; λ], parameters = z)

# @z_0: input data for z 

z₀ = [-0.08984943, -0.23977824, 0.34743529, -0.33314363, -0.01655752, -0.07218132,
-0.17829711, -0.16338844, -0.3101673, -0.25123448, 0.2592691, -0.09633576,
 0.47430885, 0.38384785, -0.15680448, 0.44287663]

res = solve(C; target_parameters = z₀)

real_solutions(res)



# characterizing max/min value
r = real_solutions(res)

θ₁, θ₂, θ₅, α₃, α₄, λ = r[1]  # change based on which solution gives a sensible result

ϵ₋₋₋₋ = θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₋₋₊ = - θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₋₋₊₋ = θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₋₊₊ = - θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₋₊₋₋ = θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₊₋₊ = - θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₋₊₊₋ = θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₋₊₊₊ = - θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₋₋₋ = - θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₋₋₊ = θ₁*α₄ + θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₊₋₊₋ = - θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₋₊₊ = θ₁*α₄ - θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₊₊₋₋ = - θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₊₋₊ = θ₁*α₄ - θ₂*α₃ - θ₁*θ₅*α₃ + θ₂*θ₅*α₄
ϵ₊₊₊₋ = - θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ - θ₂*θ₅*α₄
ϵ₊₊₊₊ = θ₁*α₄ + θ₂*α₃ + θ₁*θ₅*α₃ + θ₂*θ₅*α₄

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]


(1/2)*sum((f_inv).*square.(z₀ - ϵ))


###_____________________________________________________________________________
## Homotopy Method for Least Square Code - Case 2 option 1
###_____________________________________________________________________________

# minimize equation on slide 12
square(n) = n * n
@var θ[1:5]
@var z[1:16]
@var α[3:4]
@var f₋₋₋₋ f₋₋₋₊ f₋₋₊₋ f₋₋₊₊ f₋₊₋₋ f₋₊₋₊ f₋₊₊₋ f₋₊₊₊ f₊₋₋₋ f₊₋₋₊ f₊₋₊₋ f₊₋₊₊ f₊₊₋₋ f₊₊₋₊ f₊₊₊₋ f₊₊₊₊  θ₁ θ₂ θ₅ α₃ α₄
θ = [θ₁, θ₂, θ₅, α₃, α₄]
τ = (1/2)

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

f=[f₋₋₋₋, f₋₋₋₊, f₋₋₊₋, f₋₋₊₊, f₋₊₋₋, f₋₊₋₊, f₋₊₊₋, f₋₊₊₊, f₊₋₋₋, f₊₋₋₊, f₊₋₊₋,f₊₋₊₊, f₊₊₋₋, f₊₊₋₊, f₊₊₊₋, f₊₊₊₊]


@var ϵ₋₋₋₋ ϵ₋₋₋₊ ϵ₋₋₊₋ ϵ₋₋₊₊ ϵ₋₊₋₋ ϵ₋₊₋₊ ϵ₋₊₊₋ ϵ₋₊₊₊ ϵ₊₋₋₋ ϵ₊₋₋₊ ϵ₊₋₊₋ ϵ₊₋₊₊ ϵ₊₊₋₋ ϵ₊₊₋₊ ϵ₊₊₊₋ ϵ₊₊₊₊ α₃ α₄

ϵ₋₋₋₋ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₋₊ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₊₋ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₋₊₊ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₋ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₊ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₊₋ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₊₊₊ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₋ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₊ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₊₋ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₋₊₊ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₋ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₊ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₊₋ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₊₊₊ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]

#define constraints: τ = 1/2
g = θ₁ *θ₂ - 0.5



# Homotopy Implementation

@var Lagrangian λ L
f_inv = map(inv,f)
L = (1/2)*sum((f_inv).*square.(z- ϵ))
#Lagrangian method 
Lagrangian = L + λ*g
#Jacobian Matrix 
J = differentiate(Lagrangian,θ)
C = System([J; g], variables = [θ; λ], parameters = z)

# @z_0: input data for z 

z₀ = [-0.08984943, -0.23977824, 0.34743529, -0.33314363, -0.01655752, -0.07218132,
-0.17829711, -0.16338844, -0.3101673, -0.25123448, 0.2592691, -0.09633576,
 0.47430885, 0.38384785, -0.15680448, 0.44287663]

res = solve(C; target_parameters = z₀)

real_solutions(res)



# characterizing max/min value
r = real_solutions(res)

θ₁, θ₂, θ₅, α₃, α₄, λ = r[1]  # change based on which solution gives a sensible result

ϵ₋₋₋₋ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₋₊ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₋₊₋ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₋₊₊ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₋ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₋₊ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₋₊₊₋ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₋₊₊₊ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₋ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₋₊ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₋₊₋ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₋₊₊ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₋ = - θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₋₊ = - θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ - θ₂*θ₅*α₃
ϵ₊₊₊₋ = θ₁*θ₅*α₃ - θ₂*θ₅*α₄ - θ₁*θ₅*α₄ + θ₂*θ₅*α₃
ϵ₊₊₊₊ = θ₁*θ₅*α₃ + θ₂*θ₅*α₄ + θ₁*θ₅*α₄ + θ₂*θ₅*α₃

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]


(1/2)*sum((f_inv).*square.(z₀ - ϵ))


###_____________________________________________________________________________
## Homotopy Method for Least Square Code - Case 2 option 2
###_____________________________________________________________________________

# minimize equation on slide 12
square(n) = n * n
@var θ[1:5]
@var z[1:16]
@var α[3:5]
@var f₋₋₋₋ f₋₋₋₊ f₋₋₊₋ f₋₋₊₊ f₋₊₋₋ f₋₊₋₊ f₋₊₊₋ f₋₊₊₊ f₊₋₋₋ f₊₋₋₊ f₊₋₊₋ f₊₋₊₊ f₊₊₋₋ f₊₊₋₊ f₊₊₊₋ f₊₊₊₊  θ₁ θ₂ α₅ α₃ α₄
θ = [θ₁, θ₂, α₅, α₃, α₄]
τ = (1/2)

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

f=[f₋₋₋₋, f₋₋₋₊, f₋₋₊₋, f₋₋₊₊, f₋₊₋₋, f₋₊₋₊, f₋₊₊₋, f₋₊₊₊, f₊₋₋₋, f₊₋₋₊, f₊₋₊₋,f₊₋₊₊, f₊₊₋₋, f₊₊₋₊, f₊₊₊₋, f₊₊₊₊]


@var ϵ₋₋₋₋ ϵ₋₋₋₊ ϵ₋₋₊₋ ϵ₋₋₊₊ ϵ₋₊₋₋ ϵ₋₊₋₊ ϵ₋₊₊₋ ϵ₋₊₊₊ ϵ₊₋₋₋ ϵ₊₋₋₊ ϵ₊₋₊₋ ϵ₊₋₊₊ ϵ₊₊₋₋ ϵ₊₊₋₊ ϵ₊₊₊₋ ϵ₊₊₊₊ α₃ α₄

ϵ₋₋₋₋ = α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₋₋₊ = - α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₋₋₊₋ = - α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₋₋₊₊ = α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₊₋₋ = α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₋₊₋₊ = - α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₊₊₋ = - α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₊₊₊ = α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₋₋₋ = α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₋₋₊ = - α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₊₋₊₋ = - α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₊₋₊₊ = α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₊₋₋ = α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₊₊₋₊ = - α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₊₊₋ = - α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₊₊₊ = α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]

#define constraints: τ = 1/2
g = θ₁ *θ₂ - 0.5



# Homotopy Implementation

@var Lagrangian λ L
f_inv = map(inv,f)
L = (1/2)*sum((f_inv).*square.(z- ϵ))
#Lagrangian method 
Lagrangian = L + λ*g
#Jacobian Matrix 
J = differentiate(Lagrangian,θ)
C = System([J; g], variables = [θ; λ], parameters = z)

# @z_0: input data for z 

z₀ = [-0.08984943, -0.23977824, 0.34743529, -0.33314363, -0.01655752, -0.07218132,
-0.17829711, -0.16338844, -0.3101673, -0.25123448, 0.2592691, -0.09633576,
 0.47430885, 0.38384785, -0.15680448, 0.44287663]

res = solve(C; target_parameters = z₀)

real_solutions(res)



# characterizing max/min value
r = real_solutions(res)

θ₁, θ₂, α₅, α₃, α₄, λ = r[6]  # change based on which solution gives a sensible result

ϵ₋₋₋₋ = α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₋₋₊ = - α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₋₋₊₋ = - α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₋₋₊₊ = α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₊₋₋ = α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₋₊₋₊ = - α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₊₊₋ = - α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₋₊₊₊ = α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₋₋₋ = α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₋₋₊ = - α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₊₋₊₋ = - α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₊₋₊₊ = α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₊₋₋ = α₃*α₄ - θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ - θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃
ϵ₊₊₋₊ = - α₃*α₄ - θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ - θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₊₊₋ = - α₃*α₄ + θ₁*α₅*α₃ - θ₂*α₅*α₄ - θ₁*α₅*α₄ + θ₂*α₅*α₃ - θ₂*θ₁*α₅*α₃
ϵ₊₊₊₊ = α₃*α₄ + θ₁*α₅*α₃ + θ₂*α₅*α₄ + θ₁*α₅*α₄ + θ₂*α₅*α₃ + θ₂*θ₁*α₅*α₃

ϵ = [ϵ₋₋₋₋, ϵ₋₋₋₊, ϵ₋₋₊₋, ϵ₋₋₊₊, ϵ₋₊₋₋, ϵ₋₊₋₊, ϵ₋₊₊₋, ϵ₋₊₊₊, ϵ₊₋₋₋, ϵ₊₋₋₊, ϵ₊₋₊₋, ϵ₊₋₊₊, ϵ₊₊₋₋, ϵ₊₊₋₊, ϵ₊₊₊₋, ϵ₊₊₊₊]


(1/2)*sum((f_inv).*square.(z₀ - ϵ))



### 
###_____________________________________________________________________________
## Calculation of simple boundary case for MXM project
###_____________________________________________________________________________

# This is code I used to help calculate an analytic solution to one of the
# simplest boundary cases. Specifically the boundary case of configuration 1 (ie
# topology 13|24) with boundary constraints θ₁=1 and θ₃=3.

@var x y z[1:16] # For convenience, we write x := θ₂, y := α₄. z is the vector of data.
site_patterns=Iterators.product(fill([1;-1],4)...)|>collect # Define 2x2x2x2 array
site_patterns=map(i->reverse(B_tmp[i]),1:16) # reorder the site patterns and convert to list
reverse!(site_patterns) # reorder the list to be more natural
s=site_patterns # rename it

# Define least squares cost function

# First, define the terms i=1,...,16
X(i) = 1/((1+(1/2)*s[i][1]*s[i][2])/16)*(z[i] - s[i][2]*s[i][4]*x*y-s[i][1]*s[i][4]*y/(2x))^2

# Then add them
LS=sum([X(i) for i in 1:16])
dxLS = differentiate(LS, x)
dyLS = differentiate(LS, y)

# Setting dxLS=0 and dyLS=0, it follows that x^2*dyLS=0 and x^3*dyLS/y=0 (by
# assuming x,y>0). This multiplication eliminates stuff in the denominators, so
# we obtain two polynomial equations:
HomotopyContinuation.expand(x^2*dyLS) # quaric in x, linear in y
HomotopyContinuation.expand(x^3*dxLS/y) # quartic in x, linear in y

# The results, after simplifying and grouping terms by hand:
#
# julia> HomotopyContinuation.expand(x^2*dyLS)
# (512/3)*y
# + 32*(z₁₁ + z₆ + z₈ + z₉ - z₁₀ - z₁₂ - z₅ - z₇)*x
# + (32/3) * (z₁₃ + z₁₅ + z₂ + z₄ - z₁₄ - z₁₆ - z₁ - z₃)*x
# - (1024/3)*x^2*y 
# + 64*(z₁₀ + z₁₂ + z₅ + z₇ - z₁₁ - z₆ - z₈ - z₉)*x^3
# + (64/3)*(z₁₃ + z₁₅ + z₂ + z₄ - z₁₄ - z₁₆ - z₁ - z₃)*x^3
# + (2048/3)*x^4*y
#
# julia> HomotopyContinuation.expand(x^3*dxLS/y)
# -(512/3)*y
# + 32*(z₁₀+ z₁₂ + z₅ + z₇ - z₁₁ - z₆ - z₈ - z₉)*x
# + 32/3*(z₁ + z₁₄ + z₁₆ + z₃ - z₁₃ - z₁₅ - z₂ - z₄)*x
# + 64*(z₁₀ + z₁₂ + z₅ + z₇ - z₁₁ - z₆ - z₈ - z₉)*x^3
# + (64/3)*(z₁₃ + z₁₅ + z₂ + z₄ - z₁ - z₁₄ - z₁₆ - z₃)*x^3
# + (2048/3)*x^4*y

# Next, write these results into expressions:
@var z₁ z₂ z₃ z₄ z₅ z₆ z₇ z₈ z₉ z₁₀ z₁₁ z₁₂ z₁₃ z₁₄ z₁₅ z₁₆

expr1=(512/3)*y+32*(z₁₁+z₆+z₈+z₉-z₁₀-z₁₂-z₅-z₇)*x+(32/3)*(z₁₃+z₁₅+z₂+z₄-z₁₄-z₁₆-z₁-z₃)*x-(1024/3)*x^2*y+64*(z₁₀+z₁₂+z₅+z₇-z₁₁-z₆-z₈-z₉)*x^3+(64/3)*(z₁₃+z₁₅+z₂+z₄-z₁₄-z₁₆-z₁-z₃)*x^3+(2048/3)*x^4*y

expr2=-(512/3)*y+32*(z₁₀+z₁₂+z₅+z₇-z₁₁-z₆-z₈-z₉)*x+32/3*(z₁+z₁₄+z₁₆+z₃-z₁₃-z₁₅-z₂-z₄)*x+64*(z₁₀+z₁₂+z₅+z₇-z₁₁-z₆-z₈-z₉)*x^3+(64/3)*(z₁₃+z₁₅+z₂+z₄-z₁-z₁₄-z₁₆-z₃)*x^3+(2048/3)*x^4*y

# For fun, let's solve the system for x,y for random data:
my_system=System([expr1;expr2],variables=[x,y],parameters=z)
random_z=rand(Float64, 16)

all_solutions=solve(my_system;target_parameters=random_z,threading = true)
real_solutions(all_solutions)

# Result: for generic data, this has three solutions: (0,0), (a,b) and (-a,-b),
# where 0<a,b<1.
