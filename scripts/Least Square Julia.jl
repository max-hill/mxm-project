obs1 = [-1 -1 -1 -1 -1 -1 -1 -1 +1 +1 +1 +1 +1 +1 +1 +1; 
                -1 -1 -1 -1 +1 +1 +1 +1 -1 -1 -1 -1 +1 +1 +1 +1;
        -1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 +1 -1 -1 +1 +1;
    -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1 -1 +1]


tau = 0.5

k = 300

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
freqtest = rand(d)


using Ipopt,JuMP
model = Model(Ipopt.Optimizer)
@variable(model,x[1:16,1:5]>=0)
#set each one variable to be larger than 0
i = 1
for j in 1:16
   @constraint(model,x[j,i]<=1)
end
i = 2
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
i = 3
for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
end
i = 4
for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
end
i = 5
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
for j in 1:16
    @NLconstraint(model,x[j,1]*x[j,2]*x[j,5] == 0.5)
end
@NLobjective(model, Min, sum(q[i]*((0.5-obs1[1,i]*obs1[3,i]*x[i,1]*x[i,3]
                -obs1[2,i]*obs1[4,i]*x[i,2]*x[i,4]
                -obs1[1,i]*obs1[4,i]*x[i,1]*x[i,5]*x[i,4]
                -obs1[1,i]*obs1[4,i]*x[i,1]*x[i,5]*x[i,4])^2) for i in 1:16))
optimize!(model)
status = termination_status(model)
asdgf =objective_value(model)

function test1(a,b,c,d,e)
model = Model(Ipopt.Optimizer)
@variable(model,x[1:16,1:5]>=0)
#set each one variable to be larger than 0
i = 1
for j in 1:16
   @constraint(model,x[j,i]<=1)
end
if a == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end 
i = 2
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
if b == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end
i=3
    for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
    end
if c == 1
    for j in 1:16 
        @constraint(model,x[j,i]==0)
    end
end
i=4
    for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
    end
if d == 1
    for j in 1:16 
        @constraint(model,x[j,i]==0)
    end
end
i = 5
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
if e == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end
for j in 1:16
    @NLconstraint(model,x[j,1]*x[j,2]*x[j,5] == 0.5)
end
@NLobjective(model, Min, sum(q[i]*((0.5-obs1[1,i]*obs1[3,i]*x[i,1]*x[i,3]
                -obs1[2,i]*obs1[4,i]*x[i,2]*x[i,4]
                -obs1[1,i]*obs1[4,i]*x[i,1]*x[i,5]*x[i,4]
                -obs1[1,i]*obs1[4,i]*x[i,1]*x[i,5]*x[i,4])^2) for i in 1:16))
optimize!(model)
status = termination_status(model)
    return objective_value(model)
end

aset = test1(0,0,0,0,0)
aset

println(objective_value(model))
for i in 1:16
    for j in 1:5
        print(" ",round(value(x[i,j]),digits =2))
    end
    println()
    a = value(x[i,1])
    b = value(x[i,2])
    c =  value(x[i,5])
    println("the multiple of 125 is ", a*b*c)
end

 print(round(value(x[1,2]),digits =2))

using Ipopt,JuMP
model = Model(Ipopt.Optimizer)
@variable(model,x[1:16,1:5]>=0)
i = 1
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
i = 2
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
i = 5
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
for j in 1:16
    @NLconstraint(model,x[j,1]*x[j,2] == 0.5)
end
@NLobjective(model, Min, sum(q[i]*((0.5-obs1[1,i]*obs1[3,i]*x[i,1]*x[i,3]*x[i,5]
                -obs1[2,i]*obs1[4,i]*x[i,2]*x[i,4]*x[i,5]
                -obs1[1,i]*obs1[4,i]*x[i,1]*x[i,5]*x[i,4]
                -obs1[2,i]*obs1[3,i]*x[i,2]*x[i,5]*x[i,3])^2) for i in 1:16))
optimize!(model)

function test2(a,b,c,d,e)
model = Model(Ipopt.Optimizer)
@variable(model,x[1:16,1:5]>=0)
#set each one variable to be larger than 0
i = 1
for j in 1:16
   @constraint(model,x[j,i]<=1)
end
if a == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end 
i = 2
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
if b == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end
i=3
    for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
    end
if c == 1
    for j in 1:16 
        @constraint(model,x[j,i]==0)
    end
end
i=4
    for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
    end
if d == 1
    for j in 1:16 
        @constraint(model,x[j,i]==0)
    end
end
i = 5
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
if e == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end
for j in 1:16
    @NLconstraint(model,x[j,1]*x[j,2] == 0.5)
end
@NLobjective(model, Min, sum(q[i]*((0.5-obs1[1,i]*obs1[3,i]*x[i,1]*x[i,3]*x[i,5]
                -obs1[2,i]*obs1[4,i]*x[i,2]*x[i,4]*x[i,5]
                -obs1[1,i]*obs1[4,i]*x[i,1]*x[i,5]*x[i,4]
                -obs1[2,i]*obs1[3,i]*x[i,2]*x[i,5]*x[i,3])^2) for i in 1:16))
optimize!(model)
status = termination_status(model)
    return objective_value(model)
end

asdfas = test2(0,0,0,0,0)
asdfas

println(objective_value(model))
for i in 1:16
    for j in 1:5
        print(" ",round(value(x[i,j]),digits =2))
    end
    println()
    a = value(x[i,1])
    b = value(x[i,2])
    c =  value(x[i,5])
    println("the multiple of 125 is ", a*b*c)
end

using Ipopt,JuMP
model = Model(Ipopt.Optimizer)
@variable(model,x[1:16,1:5]>=0)
i = 1
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
i = 2
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
i = 5
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
for j in 1:16
    @NLconstraint(model,x[j,1]*x[j,2]*x[j,5] == 0.5)
end
@NLobjective(model, Min, sum(q[i]*((0.5-obs1[1,i]*obs1[4,i]*x[i,1]*x[i,4]
                -obs1[2,i]*obs1[3,i]*x[i,2]*x[i,3]
                -obs1[1,i]*obs1[3,i]*x[i,1]*x[i,5]*x[i,3]
                -obs1[2,i]*obs1[4,i]*x[i,2]*x[i,5]*x[i,4])^2) for i in 1:16))
optimize!(model)
status = termination_status(model)

function test3(a,b,c,d,e)
model = Model(Ipopt.Optimizer)
@variable(model,x[1:16,1:5]>=0)
#set each one variable to be larger than 0
i = 1
for j in 1:16
   @constraint(model,x[j,i]<=1)
end
if a == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end 
i = 2
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
if b == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end
i=3
    for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
    end
if c == 1
    for j in 1:16 
        @constraint(model,x[j,i]==0)
    end
end
i=4
    for j in 1:16
         @constraint(model,x[j,i]<=sqrt(k))
    end
if d == 1
    for j in 1:16 
        @constraint(model,x[j,i]==0)
    end
end
i = 5
for j in 1:16
    @constraint(model,x[j,i]<=1)
end
if e == 1
    for j in 1:16 
        @constraint(model,x[j,i]==1)
    end
end
for j in 1:16
    @NLconstraint(model,x[j,1]*x[j,2]*x[j,5] == 0.5)
end
@NLobjective(model, Min, sum(q[i]*((0.5-obs1[1,i]*obs1[4,i]*x[i,1]*x[i,4]
                -obs1[2,i]*obs1[3,i]*x[i,2]*x[i,3]
                -obs1[1,i]*obs1[3,i]*x[i,1]*x[i,5]*x[i,3]
                -obs1[2,i]*obs1[4,i]*x[i,2]*x[i,5]*x[i,4])^2) for i in 1:16))
optimize!(model)
status = termination_status(model)
    return objective_value(model)
end

sadfsa = test3(0,0,0,0,0)
sadfsa

println(objective_value(model))
for i in 1:16
    for j in 1:5
        print(" ",round(value(x[i,j]),digits =2))
    end
    println()
    a = value(x[i,1])
    b = value(x[i,2])
    c =  value(x[i,5])
    println("the multiple of 125 is ", a*b*c)
end

function findsmall()
    m = Array{Float64}(undef, 3, 6)
    arrayi = [0,0,0,0,0]
    m[1,1] = test1(arrayi[1],arrayi[2],arrayi[3],arrayi[4],arrayi[5])
    m[2,1] = test2(arrayi[1],arrayi[2],arrayi[3],arrayi[4],arrayi[5])
    m[3,1] = test3(arrayi[1],arrayi[2],arrayi[3],arrayi[4],arrayi[5])
    for i in 1:5
        arrayi[i] = 1
        m[1,i+1] = test1(arrayi[1],arrayi[2],arrayi[3],arrayi[4],arrayi[5])
        m[2,i+1] = test2(arrayi[1],arrayi[2],arrayi[3],arrayi[4],arrayi[5])
        m[3,i+1] = test3(arrayi[1],arrayi[2],arrayi[3],arrayi[4],arrayi[5])
        arrayi[i] =0
    end
    min = 1000
    s = findmin(m)
    return s
end
            

dsaf = findsmall()

arrayi = [0,1,0,0,0]

 m[1,1] = test1(arrayi[1],arrayi[2],arrayi[3],arrayi[4],arrayi[5])

 m = Array{Float64}(undef, 3, 6)

s = findmin(m)


