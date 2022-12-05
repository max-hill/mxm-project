 # Import necessary packages
import Pkg
Pkg.add("PDMatsExtras")
Pkg.add("Ipopt")
Pkg.add("JuMP")
using LinearAlgebra # this package comes with Julia, no need to install it
using PDMatsExtras
using Distributions
using Ipopt,JuMP,LinearAlgebra,Distributions

###
# Define Site Patterns
σ = [-1 -1 -1 -1 -1 -1 -1 -1 +1 +1 +1 +1 +1 +1 +1 +1; 
     -1 -1 -1 -1 +1 +1 +1 +1 -1 -1 -1 -1 +1 +1 +1 +1;
     -1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 +1;
     -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1]

# The number of sites (i.e. samples)
k = 100000000
τ=1/2
###

function test1(Z,τ,q)
    # This function applies the nonlinear optimization algorithm to Case 1,
    # basic case. (IE no additional constraints)xs.

    # First define the model
    model = Model(Ipopt.Optimizer)

    # Initialize the model variables – and specify that they must be nonnegative. 
    @variable(model,x1>=0)
    @variable(model,x2>=0)
    @variable(model,x3>=0)
    @variable(model,x4>=0)
    @variable(model,x5>=0)

    # Set additional constraints (these constraints apply for all cases)
    @constraint(model,x1<=1)
    @constraint(model,x2<=1)
    @constraint(model,x3<=sqrt(k))
    @constraint(model,x4<=sqrt(k))
    @constraint(model,x5<=1)
    @NLconstraint(model,x1*x2*x5 == τ)

    # Next specify the LS optimization problem. In particular, the third
    # argument is the objective function for this specific case (i.e. it is the
    # Least Squares function modified for case-specific constraints).
    @NLobjective(model, Min, (1/2*sum(((Z[i]
                                        -σ[1,i]*σ[3,i]*x1*x3
                                        -σ[2,i]*σ[4,i]*x2*x4
                                        -σ[1,i]*σ[4,i]*x1*x5*x4
                                        -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i]
                                      for i in 1:16)))
    
    # Finds the optimum of the model
    optimize!(model)

    # Store some additional information about the opimization
    status = termination_status(model)

    # Return the optimal (i.e. minimal) value of the objective function which
    # was achieved.
    return objective_value(model)
end


# Z=[ 0.04934314 -0.82975965 -0.26306304 -0.04340806  0.25569398 -0.10590171 0.15003828 -0.07734293 -0.11092292  0.19401209  0.33528518 -0.03549099 0.32687463 -0.19351553 -0.20589946  0.55405697]

# Z=[-0.07569865 -0.20756163 -0.34526203 -0.12047282  0.10219409  0.12883786  0.28248025  0.05504007  0.05976167 n-0.03856399 -0.08039081  0.31282312 -0.46831276 -0.0255255   0.16268208  0.25796903]

# Z=[-0.03387992 -0.0866895   0.15985842  0.34824363  0.10192178  0.04614786 0.01067507  0.0722576   0.22797491  0.02374262  0.21537338 -0.11314237 -0.29763172 -0.40141176 -0.56426171  0.29082169]

# Z=[ 0.21446296 -0.30243081 -0.06903343  0.01024257 -0.1291267  -0.18174438 -0.3028259   0.00650958  0.26431551 -0.0742944  -0.06907315  0.29076228 -0.00946552 -0.04134968  0.42652044 -0.03346936]

# Z=[ 0.11754145  0.24720039  0.0185607  -0.15665281 -0.23870263 -0.07454299 -0.02028725 -0.17532972  0.27826045  0.39632972  0.20278634  0.19809686 0.04671502 -0.49806674 -0.41417643  0.07226764]

###
"This function applies the nonlinear optimization algorithm to Case 1 and
    all its subcases (depending on the inputs used)"
function my_test1(Z,τ,q,mins_θ,bounds_α₃,bounds_α₄)
    # First, import function arguments
    min_θ₁, min_θ₂, min_θ₅ = mins_θ[1], mins_θ[2], mins_θ[3]
    min_α₃, max_α₃ = bounds_α₃[1], bounds_α₃[2]
    min_α₄, max_α₄ = bounds_α₄[1], bounds_α₄[2]
    # Next, define the model, along with all variables and constraints.
    model = Model(Ipopt.Optimizer)
    @variable(model, min_θ₁≤θ₁≤1)
    @variable(model, min_θ₂≤θ₂≤1)
    @variable(model, min_α₃≤α₃≤max_α₃)
    @variable(model, min_α₄≤α₄≤max_α₄)
    @variable(model, min_θ₅≤θ₅≤1)
    @NLconstraint(model,θ₁*θ₂*θ₅ == τ)
    # Define the least squares function we seek to minimize (this formula
    # applies for case 1 and all its subcases):
    objective_function=@NLexpression(model,
                                     (1/2*sum(((Z[i]
                                                -σ[1,i]*σ[3,i]*θ₁*α₃
                                                -σ[2,i]*σ[4,i]*θ₂*α₄
                                                -σ[1,i]*σ[4,i]*θ₁*θ₅*α₄
                                                -σ[2,i]*σ[3,i]*θ₂*θ₅*α₃)^2)/q[i]
                                              for i in 1:16)))
    # Now do the optimization
    @NLobjective(model, Min, objective_function)
    optimize!(model)
    # Return the result in the form [score, [coordinates], [bounds]]
    return [objective_value(model),
            [value(θ₁), value(θ₂), value(α₃), value(α₄), value(θ₅)],
            [mins_θ, bounds_α₃, bounds_α₄]]
end
###

###
# Compute the scores for all 63 subcases of case 1.
res = [my_test1(Z,1/2,q,mins_θ,bounds_α₃,bounds_α₄)
       for mins_θ in [[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0]]
           for bounds_α₃ in [[0,0], [0,sqrt(k)], [sqrt(k),sqrt(k)]]
               for bounds_α₄ in [[0,0], [0,sqrt(k)], [sqrt(k),sqrt(k)]]]

length(scores)
@show scores
scores = map(i->res[i][1], 1:length(res))
argmin(scores)
res[argmin(scores)][3]
###


function test11(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x5*x4
                -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test12(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*x4
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
    asdgf =objective_value(model)
    value(x1)
    value(x5)
    return objective_value(model)
end

function test13(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
    asdgf =objective_value(model)
    val = [value(x1),value(x2),value(x3),value(x4),value(x5)]
    return objective_value(model)
end

function test14(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test15(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test111(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x5*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test112(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test113(Z,τ,q)
    model = Model(Ipopt.Optimizer)
    @variable(model,x1>=0)
    @variable(model,x4>=0)
    @variable(model,x5>=0)
    #set each one variable to be larger than 0
    @constraint(model,x1<=1)
    @constraint(model,x4<=sqrt(k))
    @constraint(model,x5<=1)
    @NLconstraint(model,x1*x5 == τ)
    @NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                                       -σ[2,i]*σ[4,i]*x4
                                       -σ[1,i]*σ[4,i]*x1*x5*x4
                                       -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
    optimize!(model)
    status = termination_status(model)
    asdgf =objective_value(model)
    return objective_value(model)
end

function test114(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test115(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test116(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test117(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*x4
                -σ[1,i]*σ[4,i]*x5*x4
                -σ[2,i]*σ[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test118(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test119(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1111(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x4
                -σ[1,i]*σ[4,i]*x5*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1112(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1113(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1114(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1115(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1116(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end


function test2a1(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*x5*x4
                -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a2(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a3(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*0*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x2*x5*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a4(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x2*0*x5
                -σ[1,i]*σ[4,i]*x1*x5*0
                -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a5(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a11(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*x5*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a12(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*0*x5
                -σ[2,i]*σ[4,i]*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x5*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a13(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3*x5
                -σ[2,i]*σ[4,i]*x2*0*x5
                -σ[1,i]*σ[4,i]*x5*0
                -σ[2,i]*σ[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a14(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a15(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a16(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1== τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a17(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*0
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a18(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x2*x5*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x1*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b1(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*x5*x4
                -σ[2,i]*σ[3,i]*x2*x5*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b2(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x5*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x1*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b3(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*0*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*x2*x5*0
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x1*x2*0*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b4(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x2*0*x5
                -σ[1,i]*σ[4,i]*x1*x5*0
                -σ[2,i]*σ[3,i]*x2*x5*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x1*x2*x3*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b5(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*1
                -σ[2,i]*σ[4,i]*x2*x4*1
                -σ[1,i]*σ[4,i]*x1*1*x4
                -σ[2,i]*σ[3,i]*x2*1*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x1*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b11(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0) 
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4*x5
                -σ[1,i]*σ[4,i]*1*x5*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x2*0*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b12(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x4*x5
                -σ[1,i]*σ[4,i]*x1*x5*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b13(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*0
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x1*x2*0*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b14(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3*x5
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x5*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b15(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0) 
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*0*x5
                -σ[1,i]*σ[4,i]*x1*x5*0
                -σ[2,i]*σ[3,i]*x5*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b16(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0) 
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b17(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b18(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b21(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4
                -σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b22(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*0
                -σ[2,i]*σ[4,i]*x4
                -σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b23(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b24(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0) 
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0
                -σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x3
                -σ[1,i]*σ[2,i]*σ[3,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test3(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test31(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[3,i]*x3*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test32(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x3
                -σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test33(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[3,i]*x1*0*x5
                -σ[2,i]*σ[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test34(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test35(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*x2*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test311(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*x3
                -σ[1,i]*σ[3,i]*x3*x5
                -σ[2,i]*σ[4,i]*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test312(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test313(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[3,i]*x3*x5
                -σ[2,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test314(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x4
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[3,i]*x3
                -σ[2,i]*σ[4,i]*x2*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test315(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test316(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x3
                -σ[1,i]*σ[3,i]*x1*x3*x5
                -σ[2,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test317(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*x3
                -σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test318(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*x1*x4
                -σ[2,i]*σ[3,i]*0
                -σ[1,i]*σ[3,i]*0
                -σ[2,i]*σ[4,i]*x2*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test319(Z,τ,q)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == τ)
@NLobjective(model, Min, 1/2*sum(((Z[i]-σ[1,i]*σ[4,i]*0
                -σ[2,i]*σ[3,i]*x2*x3
                -σ[1,i]*σ[3,i]*x1*x3
                -σ[2,i]*σ[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function getforonetime(Z,τ,q)
    testwhole = zeros(4,25) .+100 # Matrix to store results (default value of
                                  # 100 which is so large it will be ignored).
                                  # The values of this matrix will store the
                                  # values of the objective function (The Least
                                  # Squares score – we want to find the case
                                  # with the smallest value)
    testwhole = settest1(testwhole,Z,τ,q) # Run all subcases of Case 1.
    testwhole = settest2(testwhole,Z,τ,q) # Run all subcases of Case 2, option 1
                                          
    # testwhole = settest3(testwhole,Z,τ,q) # Run all subcases Case 2, option 2
    # (Commented out for now so we don't redo the same subcases of case 2 in the
    # same run.)
    
    testwhole = settest4(testwhole,Z,τ,q) # Run all subcasese of Case 3
    return testwhole
end

function settest1(testwhole,Z,τ,q)
    # Go through each subcase of Case 1. 
    testwhole[1,1] = test1(Z,τ,q)
    testwhole[1,2] = test11(Z,τ,q)
    testwhole[1,3] = test12(Z,τ,q)
    testwhole[1,4] = test13(Z,τ,q)
    testwhole[1,5] = test14(Z,τ,q)
    testwhole[1,6] = test15(Z,τ,q)
    testwhole[1,7] = test111(Z,τ,q)
    testwhole[1,8] = test112(Z,τ,q)
    testwhole[1,9] = test113(Z,τ,q)
    testwhole[1,10] = test114(Z,τ,q)
    testwhole[1,11] = test115(Z,τ,q)
    testwhole[1,12] = test116(Z,τ,q)
    testwhole[1,13] = test117(Z,τ,q)
    testwhole[1,14] = test118(Z,τ,q)
    testwhole[1,15] = test1111(Z,τ,q)
    testwhole[1,16] = test1112(Z,τ,q)
    testwhole[1,17] = test1113(Z,τ,q)
    testwhole[1,18] = test1114(Z,τ,q)
    testwhole[1,19] = test1115(Z,τ,q)
    testwhole[1,20] = test1116(Z,τ,q)
    return testwhole
end

function settest2(testwhole,Z,τ,q)
    testwhole[2,1] = test2a(Z,τ,q)
    testwhole[2,2] = test2a1(Z,τ,q)
    testwhole[2,3] = test2a2(Z,τ,q)
    testwhole[2,4] = test2a3(Z,τ,q)
    testwhole[2,5] = test2a4(Z,τ,q)
    testwhole[2,6] = test2a5(Z,τ,q)
    testwhole[2,7] = test2a11(Z,τ,q)
    testwhole[2,8] = test2a12(Z,τ,q)
    testwhole[2,9] = test2a13(Z,τ,q)
    testwhole[2,10] = test2a14(Z,τ,q)
    testwhole[2,11] = test2a15(Z,τ,q)
    testwhole[2,12] = test2a16(Z,τ,q)
    testwhole[2,13] = test2a17(Z,τ,q)
    testwhole[2,14] = test2a18(Z,τ,q)
    return testwhole
end

function settest3(testwhole,Z,τ,q)
    testwhole[3,1] = test2b(Z,τ,q)
    testwhole[3,2] = test2b1(Z,τ,q)
    testwhole[3,3] = test2b2(Z,τ,q)
    testwhole[3,4] = test2b3(Z,τ,q)
    testwhole[3,5] = test2b4(Z,τ,q)
    testwhole[3,6] = test2b11(Z,τ,q)
    testwhole[3,7] = test2b12(Z,τ,q)
    testwhole[3,8] = test2b13(Z,τ,q)
    testwhole[3,9] = test2b14(Z,τ,q)
    testwhole[3,10] = test2b15(Z,τ,q)
    testwhole[3,11] = test2b16(Z,τ,q)
    testwhole[3,12] = test2b17(Z,τ,q)
    testwhole[3,13] = test2b18(Z,τ,q)
    testwhole[3,14] = test2b21(Z,τ,q)
    testwhole[3,16] = test2b22(Z,τ,q)
    testwhole[3,17] = test2b23(Z,τ,q)
    testwhole[3,18] = test2b24(Z,τ,q)
    return testwhole
end

function settest4(testwhole,Z,τ,q)
    testwhole[4,1] = test3(Z,τ,q)
    testwhole[4,2] = test31(Z,τ,q)
    testwhole[4,3] = test32(Z,τ,q)
    testwhole[4,4] = test33(Z,τ,q)
    testwhole[4,5] = test34(Z,τ,q)
    testwhole[4,6] = test35(Z,τ,q)
    testwhole[4,7] = test311(Z,τ,q)
    testwhole[4,8] = test312(Z,τ,q)
    testwhole[4,9] = test313(Z,τ,q)
    testwhole[4,10] = test314(Z,τ,q)
    testwhole[4,11] = test315(Z,τ,q)
    testwhole[4,12] = test316(Z,τ,q)
    testwhole[4,13] = test317(Z,τ,q)
    testwhole[4,14] = test318(Z,τ,q)
    return testwhole
end


count = zeros(4,25)



function setcount(τ)
    # The data will be generated on a quartet with leaves 1,2,3, and 4, in which
    # leaves 3 and 4 are infinitely far away and leaves 1 and 2 are distance τ=1/2
    # from each other. We call this the DATA GENERATION SETTING.

    # Compute the true probability of each site pattern under the data generation
    # setting.
    q = [ (1/16) * (1 + σ[1,j] * σ[2,j] * τ) for j in 1:16]
    
    # We will approximate date generated under the DATA GENERATION SETTING using
    # data generated by a multivariate gaussian (by CLT this is good approxmation
    # for k large). The appropriate mean and covariance matrix for this
    # approximation are defined next.
    mean = k*q
    cov = PSDMat( k*(diagm(q) - q*q'))
    d = MvNormal(mean, cov) # Our multivariate gaussian which approximates
    # large-k data
    
    # Generate a the data vector (which takes the form of the deviation of
    # sample frequencies from true probabilities)
    Z = (rand(d)-mean)/sqrt(k)
    
    testwhole = getforonetime(Z,τ,q)
    sat = findmin(testwhole)
    count[sat[2]] = count[sat[2]]+1
    return testwhole
end

for i in 1:10
    sda = setcount(0.5)
end

count


