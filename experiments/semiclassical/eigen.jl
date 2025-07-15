
"""
    Solves the eigenvalue associated with the generator -L = ∇V⋅V - β⁻¹Δ, or with 
"""
function eigenproblem(β,N,X,Ωind,U; witten=false,nev=2,verbose = false,ΔV = 0,L=nothing)
    if L === nothing
        dx = (last(X)-first(X))/(N-1)
        
        @inline prev(i) = mod(i-1,1:N)
        @inline next(i) = mod(i+1,1:N)
        
        rate(i,j) = exp(-β*(U(X[j])-U(X[i]))/2)
        rate_witten(i,j) = exp(-β*(U(X[j])-U(X[i])))

        rows = vcat(([i,i,i] for i=1:N)...)
        cols = vcat(([i,prev(i),next(i)] for i=1:N)...)

        fac = inv(β*dx^2) # scaling constant -- for consistent approximation of -L = ∇V⋅∇ - β⁻¹Δ

        if witten
            vals = vcat((fac*[rate(i,prev(i))+rate(i,next(i)),-1,-1] for i=1:N)...)
        else
            vals = vcat((fac*[rate(i,prev(i))+rate(i,next(i)),-rate(i,prev(i)),-rate(i,next(i))] for i=1:N)...)
        end

        vals *= exp(β*ΔV) # can be used so that λ1 gives the prefactor

        L = sparse(rows,cols,vals)
    end
    
    λs,us = eigs(L[Ωind,Ωind],sigma=0.0,nev=nev,which=:LR,tol=1e-12,maxiter=10000)

    λs = real.(λs)
    us = real.(us)

    return λs,us,L
end

"""
Eigenstates of the Schrödinger operator (-∂ₓ² + x²)/2 on (-inf,α) with Neumann conditions at -inf and Dirichlet at α
"""
function harm_alpha(α;N,nev,inf = 20,verbose = false)
    h =(α+inf)/N

    is = Int[]
    js = Int[]
    vals = Float64[]

    for i=1:N
        xi = -inf + (i-1)*h
        push!(is,i)
        push!(js,i)

        if i==1
            push!(vals,h*(xi^2/3+h*xi/6+h^2/30)/2+1/2h)
        else
            push!(vals,h*(xi^2/3+h^2/30)+1/h)
        end

        if i<N
            push!(is,i)
            push!(js,i+1)

            v = h*(xi^2/6+h*xi/6+h^2/20)/2 -1/2h
            push!(vals,v)
            push!(is,i+1)
            push!(js,i)
            push!(vals,v)
        end
    end

    H = sparse(is,js,vals)/h

    λs,us = eigs(H,sigma=0.0,nev=nev,which=:LR,tol=1e-12,maxiter=10000)

    return real.(λs),real.(us)
end