function interior_point_method(A, b, c; tol = 1e-5, max_time = 30, max_iter = 1000, verbose = false)
    el_time = 0.0
    st_time = 0.0
    st_time = time()
    A = sparse(A)
    (m,n) = size(A)
    e = ones(n)
    σ = 0.0
    γ = 0.99
    iter = 0
    status = "Ótimo"
    x = e; s = e; λ = ones(m)
    μ = dot(x,s)/n
    rc = A'* λ + s - c
    rb = A * x - b
    rxs = x.*s
    if verbose
        @printf("%4s  %8s  %8s  %8s  %8s  %8s  %8s\n", "iter", "||rb||", "||rc||", "||rxs||", "λ", "σ", "μ")
    end

    while norm([rc;rb;rxs]) > tol

        #Para calcular a direção afim
        b_red = [-rb;-rc + rxs./x]
        D_2 = -s./x
        D_2 = sparse(diagm(D_2))
        Q = sparse(full(([spzeros(m,m) A; A' D_2])))
        LU = lufact(Q)
        d = LU\b_red
        dλ_af = d[1:m]
        dx_af = d[m+1:n+m]
        ds_af = -(rxs + s.*dx_af) ./ x

        #Tamanho do passo de direção de escalonamento afim
        α_x = -1/min(minimum(dx_af./x),-1); α_x = min(1, α_x)
        α_s = -1/min(minimum(ds_af./s),-1); α_s = min(1, α_s)

        ##Para calcular a medida de centralização
        μ_aff = dot(x + α_x * dx_af,s + α_s * ds_af)/n
        σ = (μ_aff / μ)^3

        #Para calcular a direção corretora
        rxs = dx_af .* ds_af - (σ * μ) * e
        b_red = [-0*rb;0*rc + rxs./x]
        LU = lufact(Q)
        d = LU\b_red
        dλ_cor = d[1:m]
        dx_cor = d[m+1:n+m]
        ds_cor = -(rxs + s.*dx_cor) ./ x

        #Adicionando as componentes de escalonamento afim e correção
        dλ = dλ_af + dλ_cor
        dx = dx_af + dx_cor
        ds = ds_af + ds_cor

        #Tamanho do passo da soma das direções
        α_x = -1/min(minimum(dx./x),-1); α_x = min(1, γ * α_x)
        α_s = -1/min(minimum(ds./s),-1); α_s = min(1, γ * α_s)

        #Próximo ponto
        x += α_x * dx
        s += α_s * ds
        λ += α_s * dλ

        #Medida de desejabilidade
        μ = dot(x,s)/n

        rc = A' * λ + s - c
        rb = A * x - b
        rxs = x .* s

        #Imprimindo valores e parametros de saída
        iter = iter + 1
        verbose && @printf("%2d  %9.1e  %9.1e  %9.1e  %9.1e  %9.1e  %9.1e\n", iter, norm(rb), norm(rc), norm(rxs), norm(λ), σ, μ)
        if iter >= max_iter
            status = "Número máximo de iterações"
            break
        end
        el_time = time() - st_time
        if el_time >= max_time
            status = "Tempo máximo"
            break
        end
    end

    return x, dot(c,x), status, el_time
end
