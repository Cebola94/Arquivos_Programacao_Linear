using Clp
function clp_problema(problema; max_iter = 1000, max_time = 30)
  st_time = time()
  clp_valores = Clp.ClpCInterface
  mod = clp_valores.ClpModel()
  solve_clp = clp_valores.ClpSolve()
  
  clp_valores.read_mps(mod, problema)

  clp_valores.set_maximum_iterations(mod,max_iter)
  clp_valores.set_maximum_seconds(mod,max_time)

  clp_valores.set_presolve_type(solve_clp,1)
  clp_valores.initial_solve_with_options(mod,solve_clp)

  status = clp_valores.status(mod)
  iter = clp_valores.number_iterations(mod)
  x = clp_valores.get_col_solution(mod)
  f = clp_valores.get_obj_value(mod)
  n = clp_valores.number_cols(mod)
  m = clp_valores.number_rows(mod)

  el_time = time() - st_time
  return x ,m, n, f, status, iter, el_time
end
