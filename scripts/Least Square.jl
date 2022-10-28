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

freqtest = [ 0.05469706  
    0.36264979  
    0.0496518  
    -0.12984051 
    -0.13395214  
    0.12726297
    0.15236475  
    0.08786582  
    0.00924208 
    -0.2553283  
    -0.17056402 
    -0.19689777
    0.01828741  
    0.1594202  
    -0.21829533 
    0.08343619]

import Pkg; Pkg.add("JuMP")

using Ipopt,JuMP
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=0.999999)
@constraint(model,x[2]<=0.99999)
@constraint(model,x[3]<=sqrt(k))
@constraint(model,x[4]<=sqrt(k))
@constraint(model,x[5]<=0.999999)
@NLconstraint(model,x[1]*x[2]*x[5] == 0.5)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x[1]*x[3]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]
                -obs1[1,i]*obs1[4,i]*x[1]*x[5]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[5]*x[3])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)

function test1(a,b,c,d,e)
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=1)
if a == 1
    @constraint(model,x[1]==1)
end
if b == 1
    @constraint(model,x[2]==1)
end
if c == 1
    @constraint(model,x[3]==0)
end
if d == 1
    @constraint(model,x[4]==0)
end
if e == 1
    @constraint(model,x[5]==1)
end
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
    return objective_value(model)
end

aset = test1(0,0,0,0,1)
aset

println(objective_value(model))
for i in 1:5
    print(" ",value(x[i]))
    println()
end

using Ipopt,JuMP
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=1)
@constraint(model,x[2]<=1)
@constraint(model,x[3]<=sqrt(k))
@constraint(model,x[4]<=sqrt(k))
@constraint(model,x[5]<=1)
@NLconstraint(model,x[1]*x[2] == 0.5)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x[1]*x[3]*x[5]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]*x[5]
                -obs1[1,i]*obs1[4,i]*x[1]*x[5]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[5]*x[3])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)

function test2a(a,b,c,d,e)
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=1)
if a == 1
    @constraint(model,x[1]==1)
end
if b == 1
    @constraint(model,x[2]==1)
end
if c == 1
    @constraint(model,x[3]==0)
end
if d == 1
    @constraint(model,x[4]==0)
end
if e == 1
    @constraint(model,x[5]==1)
end
@constraint(model,x[2]<=1)
@constraint(model,x[3]<=sqrt(k))
@constraint(model,x[4]<=sqrt(k))
@constraint(model,x[5]<=1)
@NLconstraint(model,x[1]*x[2] == 0.5)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x[1]*x[3]*x[5]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]*x[5]
                -obs1[1,i]*obs1[4,i]*x[1]*x[5]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[5]*x[3])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

asdfas = test2a(0,0,0,0,0)
asdfas

using Ipopt,JuMP
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=1)
@constraint(model,x[2]<=1)
@constraint(model,x[3]<=sqrt(k))
@constraint(model,x[4]<=sqrt(k))
@constraint(model,x[5]<=1)
@NLconstraint(model,x[1]*x[2] == 0.5)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x[1]*x[3]*x[5]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]*x[5]
                -obs1[1,i]*obs1[4,i]*x[1]*x[5]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[5]*x[3]
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x[1]*x[2]*x[3]*x[4])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)

function test2b(a,b,c,d,e)
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=1)
if a == 1
    @constraint(model,x[1]==1)
end
if b == 1
    @constraint(model,x[2]==1)
end
if c == 1
    @constraint(model,x[3]==0)
end
if d == 1
    @constraint(model,x[4]==0)
end
if e == 1
    @constraint(model,x[5]==1)
end
@constraint(model,x[2]<=1)
@constraint(model,x[3]<=sqrt(k))
@constraint(model,x[4]<=sqrt(k))
@constraint(model,x[5]<=1)
@NLconstraint(model,x[1]*x[2] == 0.5)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[3,i]*x[1]*x[3]*x[5]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]*x[5]
                -obs1[1,i]*obs1[4,i]*x[1]*x[5]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[5]*x[3]
                -obs1[1,i]*obs1[2,i]*obs1[3,i]*obs1[4,i]*x[1]*x[2]*x[3]*x[4])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

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
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x[1]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[3]
                -obs1[1,i]*obs1[3,i]*x[1]*x[3]*x[5]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]*x[5])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)

function test3(a,b,c,d,e)
model = Model(Ipopt.Optimizer)
@variable(model,x[1:5]>=0)
#set each one variable to be larger than 0
@constraint(model,x[1]<=1)
if a == 1
    @constraint(model,x[1]==1)
end
if b == 1
    @constraint(model,x[2]==1)
end
if c == 1
    @constraint(model,x[3]==0)
end
if d == 1
    @constraint(model,x[4]==0)
end
if e == 1
    @constraint(model,x[5]==1)
end
@constraint(model,x[2]<=1)
@constraint(model,x[3]<=sqrt(k))
@constraint(model,x[4]<=sqrt(k))
@constraint(model,x[5]<=1)
@NLconstraint(model,x[1]*x[2]*x[5] == 0.5)
@NLobjective(model, Min, 1/2*sum(((freqtest[i]-obs1[1,i]*obs1[4,i]*x[1]*x[4]
                -obs1[2,i]*obs1[3,i]*x[2]*x[3]
                -obs1[1,i]*obs1[3,i]*x[1]*x[3]*x[5]
                -obs1[2,i]*obs1[4,i]*x[2]*x[4]*x[5])^2)/q[i] for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)
    return objective_value(model)
end

testwhole = zeros(4,6)

control = zeros(5)

testwhole[1,1] = test1(0,0,0,0,0)
testwhole[2,1] = test2a(0,0,0,0,0)
testwhole[3,1] = test2b(0,0,0,0,0)
testwhole[4,1] = test3(0,0,0,0,0)
for i = 1:5
    control[i] = 1
    testwhole[1,i+1] = test1(control[1],control[2],control[3],control[4],control[5])
    testwhole[2,i+1] = test2a(control[1],control[2],control[3],control[4],control[5])
    testwhole[3,i+1] = test2b(control[1],control[2],control[3],control[4],control[5])
    testwhole[4,i+1] = test3(control[1],control[2],control[3],control[4],control[5])
    control[i] = 0
end

    

testwhole

findmin(testwhole)

function testforall()
    testwhole[1,1] = test1(0,0,0,0,0)
    testwhole[2,1] = test2a(0,0,0,0,0)
    testwhole[3,1] = test2b(0,0,0,0,0)
    testwhole[4,1] = test3(0,0,0,0,0)
    for i = 1:5
        control[i] = 1
        testwhole[1,i+1] = test1(control[1],control[2],control[3],control[4],control[5])
        testwhole[2,i+1] = test2a(control[1],control[2],control[3],control[4],control[5])
        testwhole[3,i+1] = test2b(control[1],control[2],control[3],control[4],control[5])
        testwhole[4,i+1] = test3(control[1],control[2],control[3],control[4],control[5])
        control[i] = 0
    end
    return findmin(testwhole)
end

testforall()
