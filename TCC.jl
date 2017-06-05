# Utilizando JuMP, MathProgBase e Clp para ler o arquivo e carregar as informações necessárias
using Clp, JuMP, MathProgBase

# Inserindo arquivos onde estão o método a ser utilizado para os testes
include("pontos_interiores_problema.jl")
include("clp_problema.jl")

metodos = [pontos_interiores,clp_problema,clp_problema_sem_presolve]
metodos = [pontos_interiores]
prob_tam = length(readdir("MPS"))
met_tam = length(metodos)


P = -ones(prob_tam, met_tam)

for (i_m,met) in enumerate(metodos)

  open("testes_$met.txt", "w") do file
    str = @sprintf("%12s  %12s  %10s  %10s  %10s  %10s  %10s\n","Problema", "dot(c,x)", "Restrições", "Variáveis", "Saida", "Iterações", "Tempo")
    print(str)
    print(file,str)
    for (i_p,prob) in enumerate(readdir("MPS"))
      if met == pontos_interiores
        try
          A, b, c, m, n = ler_problema("MPS/"prob)
          if A != []
            x, m, n, f, status, iter, tempo = met(A, b, c, m, n)
            if isnan(f) || f == Inf || f == -Inf
              status = 3
            end
          else
            x, m, n, f, status, iter, tempo = Inf, m, n, Inf, Inf, Inf, Inf
          end
        catch
          error("Problema $prob inválido")
        end
      else
        x, m, n , f, status, iter, tempo = met("MPS/"prob)
      end
      if isnan(f) || f == Inf || f == -Inf
        status = 3
      end
      if status == 0 && all(x .>= 0)
        P[i_p,i_m] = tempo
      end
      str = @sprintf("%12s  %12.5e  %10d  %10d  %10d  %10d  %10.4f\n", prob, f, m, n, status, iter, tempo)
      print(str)
      print(file,str)
    end
  end
end
