function ler(problema; lim_var = 2000, lim_rest = 2000)
  mod = Model(solver=ClpSolver())
  m_internal = MathProgBase.LinearQuadraticModel(ClpSolver())
  MathProgBase.loadproblem!(m_internal, problema)
  xlb = MathProgBase.getvarLB(m_internal)
  xub = MathProgBase.getvarUB(m_internal)
  l = MathProgBase.getconstrLB(m_internal)
  m = length(l)
  n = length(xlb)
  if all(xlb .== 0) && all(xub .== Inf) && n <= lim_var && m <= lim_rest
    u = MathProgBase.getconstrUB(m_internal)
    c = MathProgBase.getobj(m_internal)
    A = MathProgBase.getconstrmatrix(m_internal)

    n = length(l)
    b = spzeros(n)
    for i = 1:n
      if l[i] == -Inf
        b[i] = u[i]
        s = spzeros(n)
        s[i] = 1
        A = [A s]
      elseif u[i] == Inf
        b[i] = l[i]
        s = spzeros(n)
        s[i] = -1
        A = [A s]
      elseif l[i] == u[i]
        b[i] = l[i]
      else
        error("Problema $problema saiu por erro nos limitantes de u e l")
        return [],[],[]
      end
    end
    println(problema)
    m,n_1 = size(A)
    n_2 = length(c)
    c = [c;spzeros(n_1 - n_2)]
    return A,b,c
  elseif m > lim_rest || n > lim_var
    error("Problema $problema saiu por erro no tamanho do problema")
    return [],[],[]
  else
    error("Problema $problema saiu por erro nos limitantes de xub e xlb")
    return [],[],[]
  end
end
