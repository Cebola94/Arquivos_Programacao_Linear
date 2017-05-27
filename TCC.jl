# Utilizando JuMP, MathProgBase e Clp para ler o arquivo e carregar as informações necessárias
using JuMP, MathProgBase, Clp

# Utilizando BenchmarkProfiles para realizar os gráficos de performance profile
using BenchmarkProfiles, Plots

# Inserindo arquivo onde está o método a ser utilizado para os testes
include("pontos_interiores.jl")
include("clp_solve.jl")
#include("Simplex_LUfact.jl")

#=lista=["adlittle.mps","afiro.mps","agg.mps","bandm.mps","beaconfd.mps","blend.mps","brandy.mps"
,"e226.mps","israel.mps","kb2.mps","lotfi.mps",
"sc105.mps","sc205.mps","sc50a.mps","sc50b.mps","scagr25.mps","scagr7.mps","scfxm1.mps","scorpion.mps",
"sctap1.mps","share1b.mps","share2b.mps","stocfor1.mps"]
=#
#lista = ["afiro.mps","25fv47.mps","finnis.mps","dfl001.mps","qap15.mps"]
#
lista =["scsd6.mps","scsd8.mps","sctap2.mps","sctap3.mps","wood1p.mps","sctap2.mps","sctap3.mps",
"share2b.mps","kleemin3.mps","kleemin4.mps","kleemin5.mps","kleemin6.mps","kleemin7.mps","farm.mps",
"stocfor2.mps","degen3.mps","afiro.mps","blend.mps","israel.mps","sc105.mps","sc50a.mps","sc50b.mps",
"share2b.mps","stocfor1.mps","capri.mps","e226.mps","sc205.mps","sctap1.mps","stair.mps"]

metodos = [clp,interior_point_method]

prob_tam = length(lista)
met_tam = length(metodos)
println(prob_tam)
P = -ones(prob_tam, met_tam)

for (i_m,met) in enumerate(metodos)

    open("pontos_interiores$met.txt", "w") do file
        str = @sprintf("%12s  %12s  %5s  %5s\n","Problema", "dot(c,x)", "Saida", "Tempo")
        print(str)
        print(file,str)
        for (i_p,prob) in enumerate(lista)
            mod = Model(solver=ClpSolver())
            m_internal = MathProgBase.LinearQuadraticModel(ClpSolver())
            MathProgBase.loadproblem!(m_internal, "MPS/"prob)
            c = MathProgBase.getobj(m_internal)
            A = MathProgBase.getconstrmatrix(m_internal)
            m, n = size(A)
            xlb = MathProgBase.getvarLB(m_internal)
            xub = MathProgBase.getvarUB(m_internal)
            l = MathProgBase.getconstrLB(m_internal)
            u = MathProgBase.getconstrUB(m_internal)
            vtypes = MathProgBase.getvartype(m_internal)

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
                    error("Problema $prob inválido")
                end
            end
            m,n_1 = size(A)
            n_2 = length(c)
            c = [c;spzeros(n_1 - n_2)]
            try
                x ,f, status, tempo = met(A, b, c)
                if isnan(f) || f == Inf || f == -Inf
					status = 3
				end
                if status == 0
                    P[i_p,i_m] = tempo
                end
                str = @sprintf("%12s  %12.5e  %5d  %5.4f\n", prob, f, status, tempo)
                print(str)
                print(file,str)
            catch
                status = 3
                str = @sprintf("%12s  %12.5e  %5d  %5.4f\n", prob, " ", status, " ")
                print(str)
                print(file,str)
            end
        end
    end
end
performance_profile(P, ["ClpSolver", "PontosInteriores"])
png("perprof_t.png")
