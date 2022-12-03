 # Import necessary packages
Pkg.add("PDMatsExtras")
import Pkg; Pkg.add("Ipopt")
import Pkg; Pkg.add("JuMP")
using LinearAlgebra # this package comes with Julia, no need to install it
import Pkg
using PDMatsExtras
using Distributions
using Ipopt,JuMP,LinearAlgebra,Distributions


# Define Site Patterns
obs1 = [-1 -1 -1 -1 -1 -1 -1 -1 +1 +1 +1 +1 +1 +1 +1 +1; 
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
    @NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                                       -obs1[2,i]*obs1[4,i]*x2*x4
                                       -obs1[1,i]*obs1[4,i]*x1*x5*x4
                                       -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
    # Finds the optimum of the model
    optimize!(model)

    # Store some additional information about the opimization
    status = termination_status(model)

    # Return the optimal (i.e. minimal) value of the objective function which
    # was achieved.
    return objective_value(model)
end

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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*0*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*0*x5
                -obs1[1,i]*obs1[4,i]*x5*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*x3*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x3*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*0*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*0*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*x3*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*1
                -obs1[2,i]*obs1[4,i]*x2*x4*1
                -obs1[1,i]*obs1[4,i]*x1*1*x4
                -obs1[2,i]*obs1[3,i]*x2*1*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*x3*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*1*x5*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x2*0*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*0*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*0
                -obs1[2,i]*obs1[3,i]*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4*x5)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4)^2)/q[i] for i in 1:16))
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
@NLobjective(model, Min, 1/2*sum(((Z[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function getforonetime(Z,τ,q)
    testwhole = zeros(4,25) .+100 # Matrix to store results (default value of
                                  # 100 which is so large it will be ignored).
                                  # The values of this matrix will be the value
                                  # of the objective function (The Least Squares
                                  # score – we want to find the case with the
                                  # smallest value)
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

###
count = zeros(4,25)



function setcount(τ)
    # The data will be generated on a quartet with leaves 1,2,3, and 4, in which
    # leaves 3 and 4 are infinitely far away and leaves 1 and 2 are distance τ=1/2
    # from each other. We call this the DATA GENERATION SETTING.

    # Compute the true probability of each site pattern under the data generation
    # setting.
    q = zeros(16)
    for j = 1:16;  q[j] = (1 / 16) * (1 + obs1[1, j] * obs1[2, j] * τ); end

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


