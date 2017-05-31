# Utilizando JuMP, MathProgBase e Clp para ler o arquivo e carregar as informações necessárias
using JuMP, MathProgBase, Clp

# Inserindo arquivos onde estão o método a ser utilizado para os testes
include("pontos_interiores_problema.jl")
include("clp_problema.jl")

#Lista de problemas a serem resolvidos
lista =["noswot.mps","scsd6.mps","scsd8.mps","wood1p.mps","sctap2.mps","sctap3.mps",
"kleemin3.mps","kleemin4.mps","kleemin5.mps","kleemin6.mps","kleemin7.mps","farm.mps",
"stocfor2.mps","afiro.mps","blend.mps","israel.mps","sc105.mps","sc50a.mps","sc50b.mps",
"share2b.mps","stocfor1.mps","sc205.mps","sctap1.mps"]

metodos = [clp_problema,pontos_interiores_problema]

prob_tam = length(lista)
met_tam = length(metodos)
P = -ones(prob_tam, met_tam)

for (i_m,met) in enumerate(metodos)

    open("pontos_interiores_$met.txt", "w") do file
        str = @sprintf("%12s  %12s  %5s  %9s  %5s\n","Problema", "dot(c,x)", "Saida", "Iterações", "Tempo")
        print(str)
        print(file,str)
        for (i_p,prob) in enumerate(lista)
            try
                x ,f, status, iter, tempo = met("MPS/"prob)
                if isnan(f) || f == Inf || f == -Inf
					status = 3
				end
                if status == 0 && all(x .>= 0)
                    P[i_p,i_m] = tempo
                end
                str = @sprintf("%12s  %12.5e  %5d  %5d  %5.4f\n", prob, f, status, iter, tempo)
                print(str)
                print(file,str)
            catch
                status = 3
                str = @sprintf("%12s  %12.5e  %5d  %5d  %5.4f\n", prob, " ", status, " ", " ")
                print(str)
                print(file,str)
            end
        end
    end
end
