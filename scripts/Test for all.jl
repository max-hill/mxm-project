obs1 = [-1 -1 -1 -1 -1 -1 -1 -1 +1 +1 +1 +1 +1 +1 +1 +1; 
                -1 -1 -1 -1 +1 +1 +1 +1 -1 -1 -1 -1 +1 +1 +1 +1;
        -1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 +1;
    -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1]


tau = 0.5

k = 1000000000


q = zeros(16)
for j = 1:16
   q[j] = (1 / 16) * (1 + obs1[1, j] * obs1[2, j] * tau)
end

q


mean = k *q

using LinearAlgebra # this package comes with Julia, no need to install it

meandia = diagm(q)

b = transpose(q)

c = q * b

d = meandia -c 

cov = k * d
# find the Covariance of the distribution

using Distributions
d = MvNormal(mean, cov)
test = rand(d)
# the basic generated data

import Pkg; Pkg.add("Ipopt")

import Pkg; Pkg.add("JuMP")

using Ipopt,JuMP
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=1)
@constraint(model,x[2]<=1)
@constraint(model,x[3]<=sqrt(k))
@constraint(model,x[4]<=sqrt(k))
@constraint(model,x[5]<=1)
@NLconstraint(model,x[1]*x[2]*x[5] == 0.5)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x[1]*x[3]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]
                -obs1[1,i]*obs1[4,i]*x[1]*x[5]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[5]*x[3])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)

function test1(freqtest,e)
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
@NLconstraint(model,x1*x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test11(freqtest,e)
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
@NLconstraint(model,x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test12(freqtest,e)
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
@NLconstraint(model,x1*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test13(freqtest,e)
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
@NLconstraint(model,x1*x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test14(freqtest,e)
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
@NLconstraint(model,x1*x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test15(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test111(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test112(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test113(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test114(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test115(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test116(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test117(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test118(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test119(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1111(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1112(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1113(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1114(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1115(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test1116(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end


function test2a1(freqtest,e)
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
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a2(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a3(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a4(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*0*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a5(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a11(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a12(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a13(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*0*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a14(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x5*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a15(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a16(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1== e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a17(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2a18(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b1(freqtest,e)
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
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b2(freqtest,e)
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
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b3(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*x2*x5*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*0*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b4(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*0*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*x3*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b5(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*1
                -obs1[2,i]*obs1[4,i]*x2*x4*1
                -obs1[1,i]*obs1[4,i]*x1*1*x4
                -obs1[2,i]*obs1[3,i]*x2*1*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b11(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0) 
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4*x5
                -obs1[1,i]*obs1[4,i]*1*x5*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*0*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b12(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b13(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x1*x2*0*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b14(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b15(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0) 
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0*x5
                -obs1[1,i]*obs1[4,i]*x1*x5*0
                -obs1[2,i]*obs1[3,i]*x5*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b16(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0) 
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b17(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b18(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x2*x3*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b21(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4
                -obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b22(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*0
                -obs1[2,i]*obs1[4,i]*x4
                -obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b23(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0) 
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test2b24(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0) 
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0
                -obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test3(freqtest,e)
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
@NLconstraint(model,x1*x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test31(freqtest,e)
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
@NLconstraint(model,x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test32(freqtest,e)
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
@NLconstraint(model,x1*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test33(freqtest,e)
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
@NLconstraint(model,x1*x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*x1*0*x5
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test34(freqtest,e)
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
@NLconstraint(model,x1*x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test35(freqtest,e)
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
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x2*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test311(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x3>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test312(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test313(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x2*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x3*x5
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test314(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x2>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x4
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x3
                -obs1[2,i]*obs1[4,i]*x2*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test315(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x4>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x4<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x4*x5)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test316(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x5>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x5<=1)
@NLconstraint(model,x1*x5 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x1*x3*x5
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test317(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x3>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x3<=sqrt(k))
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*x3
                -obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test318(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x4>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x4<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x1*x4
                -obs1[2,i]*obs1[3,i]*0
                -obs1[1,i]*obs1[3,i]*0
                -obs1[2,i]*obs1[4,i]*x2*x4)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function test319(freqtest,e)
model = Model(Ipopt.Optimizer)
@variable(model,x1>=0)
@variable(model,x2>=0)
@variable(model,x3>=0)
#set each one variable to be larger than 0
@constraint(model,x1<=1)
@constraint(model,x2<=1)
@constraint(model,x3<=sqrt(k))
@NLconstraint(model,x1*x2 == e)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*0
                -obs1[2,i]*obs1[3,i]*x2*x3
                -obs1[1,i]*obs1[3,i]*x1*x3
                -obs1[2,i]*obs1[4,i]*0)^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

function getforonetime(freqtest,e)
    testwhole = zeros(4,25)
    testwhole = settest1(testwhole,freqtest,e)
    testwhole = settest2(testwhole,freqtest,e)
    testwhole = settest3(testwhole,freqtest,e)
    testwhole = settest4(testwhole,freqtest,e)
    for i in 1:4
        for j in 1:25
            if testwhole[i,j] == 0
                testwhole[i,j] = 100
            end
        end
    end
    return testwhole
end

function settest1(testwhole,freqtest,e)
    testwhole[1,1] = test1(freqtest,e)
    testwhole[1,2] = test11(freqtest,e)
    testwhole[1,3] = test12(freqtest,e)
    testwhole[1,4] = test13(freqtest,e)
    testwhole[1,5] = test14(freqtest,e)
    testwhole[1,6] = test15(freqtest,e)
    testwhole[1,7] = test111(freqtest,e)
    testwhole[1,8] = test112(freqtest,e)
    testwhole[1,9] = test113(freqtest,e)
    testwhole[1,10] = test114(freqtest,e)
    testwhole[1,11] = test115(freqtest,e)
    testwhole[1,12] = test116(freqtest,e)
    testwhole[1,13] = test117(freqtest,e)
    testwhole[1,14] = test118(freqtest,e)
    testwhole[1,15] = test1111(freqtest,e)
    testwhole[1,16] = test1112(freqtest,e)
    testwhole[1,17] = test1113(freqtest,e)
    testwhole[1,18] = test1114(freqtest,e)
    testwhole[1,19] = test1115(freqtest,e)
    testwhole[1,20] = test1116(freqtest,e)
    return testwhole
end

function settest2(testwhole,freqtest,e)
    testwhole[2,1] = test2a(freqtest,e)
    testwhole[2,2] = test2a1(freqtest,e)
    testwhole[2,3] = test2a2(freqtest,e)
    testwhole[2,4] = test2a3(freqtest,e)
    testwhole[2,5] = test2a4(freqtest,e)
    testwhole[2,6] = test2a5(freqtest,e)
    testwhole[2,7] = test2a11(freqtest,e)
    testwhole[2,8] = test2a12(freqtest,e)
    testwhole[2,9] = test2a13(freqtest,e)
    testwhole[2,10] = test2a14(freqtest,e)
    testwhole[2,11] = test2a15(freqtest,e)
    testwhole[2,12] = test2a16(freqtest,e)
    testwhole[2,13] = test2a17(freqtest,e)
    testwhole[2,14] = test2a18(freqtest,e)
    return testwhole
end

function settest3(testwhole,freqtest,e)
    testwhole[3,1] = test2b(freqtest,e)
    testwhole[3,2] = test2b1(freqtest,e)
    testwhole[3,3] = test2b2(freqtest,e)
    testwhole[3,4] = test2b3(freqtest,e)
    testwhole[3,5] = test2b4(freqtest,e)
    testwhole[3,6] = test2b11(freqtest,e)
    testwhole[3,7] = test2b12(freqtest,e)
    testwhole[3,8] = test2b13(freqtest,e)
    testwhole[3,9] = test2b14(freqtest,e)
    testwhole[3,10] = test2b15(freqtest,e)
    testwhole[3,11] = test2b16(freqtest,e)
    testwhole[3,12] = test2b17(freqtest,e)
    testwhole[3,13] = test2b18(freqtest,e)
    testwhole[3,14] = test2b21(freqtest,e)
    testwhole[3,16] = test2b22(freqtest,e)
    testwhole[3,17] = test2b23(freqtest,e)
    testwhole[3,18] = test2b24(freqtest,e)
    return testwhole
end

function settest4(testwhole,freqtest,e)
    testwhole[4,1] = test3(freqtest,e)
    testwhole[4,2] = test31(freqtest,e)
    testwhole[4,3] = test32(freqtest,e)
    testwhole[4,4] = test33(freqtest,e)
    testwhole[4,5] = test34(freqtest,e)
    testwhole[4,6] = test35(freqtest,e)
    testwhole[4,7] = test311(freqtest,e)
    testwhole[4,8] = test312(freqtest,e)
    testwhole[4,9] = test313(freqtest,e)
    testwhole[4,10] = test314(freqtest,e)
    testwhole[4,11] = test315(freqtest,e)
    testwhole[4,12] = test316(freqtest,e)
    testwhole[4,13] = test317(freqtest,e)
    testwhole[4,14] = test318(freqtest,e)
    return testwhole
end

count = zeros(4,25)

function setcount(e)
    test123 = rand(d)
    freqtest = (test123-mean)/sqrt(k)
    testwhole = getforonetime(freqtest,e)
    sat = findmin(testwhole)
    count[sat[2]] = count[sat[2]]+1
    return testwhole
end

for i in 1:10
    setcount(0.5)
end

count


