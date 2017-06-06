# Utilizando JuMP, MathProgBase e Clp para ler o arquivo e carregar as informações necessárias
using Clp, JuMP, MathProgBase

# Inserindo arquivos onde estão o método a ser utilizado para os testes
include("pontos_interiores_problema.jl")
include("clp_problema.jl")

metodos = [pontos_interiores,clp_problema,clp_problema_sem_presolve]

lista = ["afiro.mps","agg2.mps","agg3.mps","agg.mps","bandm.mps","blend.mps",
"degen2.mps","degen3.mps","e226.mps","farm.mps",
"israel.mps","kleemin3.mps","kleemin4.mps","kleemin5.mps","kleemin6.mps","kleemin7.mps",
"kleemin8.mps","lotfi.mps","sc105.mps","sc205.mps","sc50a.mps","sc50b.mps",
"scagr25.mps","scagr7.mps","scfxm1.mps","scfxm2.mps","scfxm3.mps","scsd1.mps",
"scsd6.mps","sctap1.mps","sctap2.mps","stocfor1.mps"]

prob_tam = length(lista)
met_tam = length(metodos)

P = -ones(prob_tam, met_tam)

for (i_m,met) in enumerate(metodos)
  open("pontos_interiores_$met.txt", "w") do file
    str = @sprintf("%12s  %12s  %10s  %10s  %10s  %10s  %8s\n", "Problema", "dot(c,x)", "Restrições", "Variáveis", "Saida", "Iterações", "Tempo")
    print(str)
    print(file,str)
    for (i_p,prob) in enumerate(lista)
      if met == pontos_interiores
        A, b, c, m, n = ler_problema("MPS/"prob)
        if A != []
          x, m, n, f, status, iter, tempo = met(A, b, c, m, n)
          if isnan(f) || f == Inf || f == -Inf
            status = 3
          end
        else
          x, m, n, f, status, iter, tempo = Inf, m, n, Inf, 5, Inf, Inf
        end
      else
        x, m, n, f, status, iter, tempo = met("MPS/"prob)
      end
        str = @sprintf("%12s  %12.5e  %10d  %10d  %10d  %10d  %10.4f\n", prob, f, m, n, status, iter, tempo)
        print(str)
        print(file,str)
    end
  end
end
