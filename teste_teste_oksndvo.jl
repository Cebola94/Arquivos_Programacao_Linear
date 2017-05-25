using ForwardDiff
include("teste_oksndvo.jl")
f(x) = sin(x)
A = [1 2;2 3]
b = [1;6]
x0 = [1.0;1.0]
#c = [0;2]

x, λ, iter, ef, el_time = pre_corr_QP(f,x0,A,b)
println("x = $x, λ = $λ, iter = $iter, ef=$ef, el_time=$el_time")
