using JuMP, Clp
function clp(A, b, c)
    st_time = time()
    mod = Model(solver=ClpSolver())
    n = length(c)
    @variable(mod, x[1:n] >= 0)
    @constraint(mod, A*x .== b)
    @objective(mod, Min, dot(c,x))

    status = solve(mod)

    if string(:Optimal) == "Optimal"
        stat = 0
    else
        stat = 1
    end

    x = getvalue(x)
    f = getobjectivevalue(mod)
    el_time = time() - st_time
    return x, f, stat, el_time
end
