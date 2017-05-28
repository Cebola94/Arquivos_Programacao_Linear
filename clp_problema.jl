using JuMP, Clp
function clp_problema(problema)
    st_time = time()
    clp_valores = Clp.ClpCInterface
    mod = clp_valores.ClpModel()
    clp_valores.read_mps(mod,problema,false ,true)
    clp_valores.initial_solve(mod)
    status = clp_valores.status(mod)
    iter = clp_valores.number_iterations(mod)
    x = clp_valores.get_col_solution(mod)
    f = clp_valores.get_obj_value(mod)

    el_time = time() - st_time
    return x, f, status, iter, el_time
end
