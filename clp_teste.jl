using JuMP, Clp
function clp_solve(A, b, c)
    st_time = time()
    mod = Model(solver=ClpSolver())
    n = length(c)
    @variable(mod, x[1:n] >= 0)
    @constraint(mod, A*x .== b)
    @objective(mod, Min, dot(c,x))

    status = solve(mod)
    el_time = time() - st_time
    #iter = get_number_iteration(mod)
    x = getvalue(x)
    f = getobjectivevalue(mod)
    return x, f, status, el_time
end
