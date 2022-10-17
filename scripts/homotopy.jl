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


# system of equations

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
n = 4 # number of leaves (4 in our case)
H = differentiate(differentiate(Lagrangian,f),f) # The Hessian of the Lagrangian.

random_u = n₀ # the data

function classify_critical_point(critical_point::Array)

    "This function implements the bordered determinental criterion (i.e. 2nd
    derivative test) to test if a given real critical point is a local maximum.
    Specifically, we follow Theorem 12 on page 155 of Magnus and Neudecker, Jan
    R. and Heinz. Matrix Differential Calculus with Applications in Statistics
    and Econometrics. Second (revised) edition, 1999. It returns 1 if the point
    is a maximum, -1 if the point is a minimum, and 0 if the test is
    indeterminate."

    c=critical_point[1:2^n] # The pᵢ's
    d=critical_point[(2^n+1):2^n+m] # The λᵢ's
    M=fill(0,m,m) # Square m×m matrix of zeros.

    # First, we need to test that the the m×m upper-left square submatrix of J
    # has nonzero determinant when evaluated at the critical point. Magnus et al
    # claim this condition is necessary but that we can always relabel our
    # variables so that this condition holds. I have not implemented relabeling
    # but will if it becomes necessary to do so. For now, the function just
    # tests to see if the condition holds and gives up if it does not.

    Jₘ=view(J,:,1:m)
    Jₘ=subs(Jₘ,[p;λ;u]=>[c;d;random_u]) # evaluate it at the critical point
    Jₘ=convert(Array{Float32,2},Jₘ) # convert it to float
    if det(Jₘ)==0
        return "Error: unmet condition. Indeterminate."
    else
        tmp_vector=[]
        for r in (m+1):2^n
            Hᵣ=view(H,1:r,1:r) # r×r submatrix of the Hessian of the Lagrangian
            Jᵣ=view(J,:,1:r)   # First r columns of the Jacobian of constraints

            # Next, we construct Δ, which is a symmetric (m+r)×(m+r) minor of
            # the (m+n)×(m+n) bordered Hessian. (In particular, Δ is the matrix
            # defined by equation (6) in Magnus et al. on page 155.) We then
            # evaluate Δ at the critical point, and with the randomly-chosen
            # parameters used when constructing our system. For computational
            # reasons we then tell Julia to treat the entries of Δ as numbers,
            # rather than symbolic expressions, and compute its determinant.
            
            Δ_matrix = [M Jᵣ;Jᵣ' Hᵣ]
            Δ_matrix = subs(Δ_matrix, [p;λ;u]=>[c;d;random_u])
            Δ_matrix = convert(Array{Float32,2},Δ_matrix)
            Δ_det = prod(eigvals(Δ_matrix))
            
            # Collect the determinant into a temporary vector.
            tmp_vector=[tmp_vector;Δ_det]
        end
        
        # Finally, check the determinental condition for classification.
        if [tmp_vector[i]*(-1)^(i+m)>0 for i in 1:(2^n-m)] == fill(true, 2^n-m)
            return 1 # yes, maximum. note i=r-m
        elseif [tmp_vector[i]*(-1)^m>0 for i in 1:(2^n-m)] == fill(true, 2^n-m)
            return -1 # no, minimum
        else
            return 0 # saddle-point or degenerate extrema
        end
    end
end
