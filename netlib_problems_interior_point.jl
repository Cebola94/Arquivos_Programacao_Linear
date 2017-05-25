include("interior_point.jl")
using JuMP
using MathProgBase
using Clp

mod = Model(solver=ClpSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0, DisplayInterval=1))
m_internal = MathProgBase.LinearQuadraticModel(ClpSolver())

MathProgBase.loadproblem!(m_internal, "MPS/noswot.mps")

c = MathProgBase.getobj(m_internal)
A = MathProgBase.getconstrmatrix(m_internal)
m, n = size(A)
xlb = MathProgBase.getvarLB(m_internal)
xub = MathProgBase.getvarUB(m_internal)
l = MathProgBase.getconstrLB(m_internal)
u = MathProgBase.getconstrUB(m_internal)
vtypes = MathProgBase.getvartype(m_internal)
