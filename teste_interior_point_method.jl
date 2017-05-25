using FactCheck, Krylov
include("interior_point.jl")
include("passo_newton.jl")
include("tamanho_passo.jl")
verbose = false
facts("\nmin x1 + 2x2 s.a x1 + x2 = 1, x1,x2 >= 0") do
    A = [1 2]
    b = 1
    c = [1;1]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin x1 + 2x2 s.a 2x1 + x2 = 0.5, x1,x2 >= 0") do
    A = [2 1]
    b = 0.5
    c = [1;2]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin x1 + x2 s.a x1 + 2x2 = 3, 3x1 + x2 = 3, x1,x2 >= 0") do
    A = [1 2;3 1]
    b = [3,3]
    c = [1;1]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin x1 - x2 s.a 2x1 + x2 = 4;x1 + 3x2 = 2, x1,x2 >= 0\n") do
    A = [2 1;1 3]
    b = [4;2]
    c = [1;-1]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose, tol = 1e-10)
    println("x^{*} = $x, s^{*} = $s, λ^{*} = $λ")
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin x1 + 2x2 s.a x1 + x2 = 3;x1 + 2x2 = 2, x1,x2 >= 0") do
    A = [2 1;1 2]
    b = [3;2]
    c = [1;2]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin 2x1 + 3x2 s.a x1 + 2x2 = 4;x1 + x2 = 3, x1,x2 >= 0") do
    A = [1 2;1 1]
    b = [4;3]
    c = [2;3]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin -x1 - 2x2 s.a x1 + x2 = 1, x1,x2 >= 0") do
    A = [1 1]
    b = 1
    c = [-1;-2]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin 2x1 + 3x2 s.a x1 + x2 - x3 = 4;  x1 + x2 - x4 = 3,x1,x2,x3,x4 >= 0") do
    A = [1 1 -1 0;1 1 0 -1]
    b = [4,3]
    c = [2;2;0;0]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end

facts("\nmin 2x1 + x2 s.a x1 + x2 + x3 = 5;  x1 + x3 = 5,x1,x2,x3 >= 0") do
    A = [1 1 1;1 0 1]
    b = [5,5]
    c = [2;1;0]
    x, λ, s, cx, cx_dual, iter, ef, el_time = interior_point_method(A,b,c, verbose = verbose)
    @fact cx --> roughly(0, 1e-5)
    @fact cx_dual --> roughly(0, 1e-5)
    @fact ef --> 0
    @fact iter --> less_than(1000)
end
