# Utilizando JuMP, MathProgBase e Clp para ler o arquivo e carregar as informações necessárias
using JuMP, MathProgBase, Clp

# Inserindo arquivo onde está o método a ser utilizado para os testes
include("interior_point.jl")
include("clp_teste.jl")

lista=["25fv47.mps","80bau3b.mps","adlittle.mps","afiro.mps","agg.mps","agg2.mps","agg3.mps","bandm.mps","beaconfd.mps","blend.mps","bnl1.mps","bnl2.mps","boeing1.mps","boeing2.mps","bore3d.mps","brandy.mps","capri.mps","cycle.mps","czprob.mps","d2q06c.mps","d6cube.mps","degen2.mps","degen3.mps","dfl001.mps"
,"e226.mps","etamacro.mps","fffff800.mps","finnis.mps","fit1d.mps","fit1p.mps","fit2d.mps","fit2p.mps","forplan.mps","ganges.mps","gfrd-pnc.mps","greenbea.mps","greenbeb.mps","grow15.mps","g22.mps","grow7.mps","israel.mps","kb2.mps","lotfi.mps","maros.mps","maros-r7.mps","modszk1.mps","nesm.mps","perold.mps","pilot.mps",
"pilot.ja.mps","pilot.we.mps","pilot4.mps","pilot87.mps","pilotnov.mps","qap8.mps","qap12.mps","qap15.mps","recipe.mps","sc105.mps","sc205.mps","sc50a.mps","sc50b.mps","scagr25.mps","scagr7.mps","scfxm1.mps","scfxm2.mps","scfxm3.mps","scorpion.mps","scrs8.mps","scsd1.mps","scsd6.mps","scsd8.mps","sctap1.mps","sctap2.mps",
"sctap3.mps","seba.mps","share1b.mps","share2b.mps","shell.mps","ship04l.mps","ship04s.mps","ship08l.mps","ship08s.mps","ship12l.mps","ship12s.mps","sierra.mps","stair.mps","standata.mps","standgub.mps","standmps.mps","stocfor1.mps","stocfor2.mps","stocfor3.mps","truss.mps","tuff.mps","vtp.base.mps","wood1p.mps","woodw.mps"]

lista = ["afiro.mps","25fv47.mps"]

#open("log_penalidade_caixa.txt", "w") do file
    #str = @sprintf("%12s  %10s  %10s  %10s\n","Problema", "dot(c,x)", "dot(λ,b)","rb","rc","iter","ef","el_time")
    #print(str)
    #print(file,str)

    for problema in lista

        mod = Model(solver=ClpSolver())
        m_internal = MathProgBase.LinearQuadraticModel(ClpSolver())
        MathProgBase.loadproblem!(m_internal, problema)
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
                error("Problema inválido")
            end
            if 
        end
        m,n_1 = size(A)
        n_2 = length(c)
        c = [c;spzeros(n_1 - n_2)]
        x ,f_clp, status_clp, eltime_clp = clp_solve(A, b, c)
        x ,f_pi, status_pi, eltime_pi = interior_point_method(A, b, c)
        println("Valor ótimo Clp= $f_clp, Status Clp = $status_clp, Tempo Clp = $eltime_clp")
        println("Valor ótimo PI= $f_pi, Status PI = $status_pi, Tempo Clp = $eltime_pi")
        #str = @sprintf("%12s  %10.4e  %10.4e  %10.4e\n",problema, f_pi, iter, el_time)
        #print(str)
        #print(file,str)
    end
#end
