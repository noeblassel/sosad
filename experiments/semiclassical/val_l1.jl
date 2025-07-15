using SparseArrays,Arpack,LinearAlgebra,Zygote,Distributions,ProgressMeter,LaTeXStrings,Plots
include("eigen.jl")
include("pot_opt.jl")
# V(x) = 1-cos((x-0.1)/(1+0.5x^2)) + 0.06507069774709273*(x-1.3177)/2.8354


ab_pairs = [(0.5,0.3),(1.0,-0.3),(0.0,0.0)]
markershapes=[:circle,:xcross,:cross]
styles=[:dashdot,:dot,:dash]
colors = [:green,:blue,:red]
l1_plot = plot(ylabel=L"\lambda_1",xlabel=L"\beta")
βrange = range(0,10,21)

for (i,t)=enumerate(ab_pairs)
    a,b = t
    λ1s = Float64[]
    
    @showprogress for β = βrange
        Ωind = Ω_f(β,a,b)
        try
            λs, us = eigenproblem(β,N,X,Ωind,V;nev=2,witten=true)
            push!(λ1s,first(λs))
        catch
            push!(λ1s,NaN)
        end
    end
    
    println(λ1s)
    scatter!(l1_plot,βrange,λ1s,label=L"\alpha^{(1)}=%$a,\alpha^{(2)}=%$b",markershape=markershapes[i],color=colors[i],markerstrokecolor=colors[i],markercolor=((colors[i] == :green) ? :white : colors[i]))
    Φ(x) = cdf(Normal(),Float64(x))
    prefactor = (sqrt(ν1)/Φ(sqrt(ν1)*a) + sqrt(ν2)/Φ(sqrt(ν2)*b))*sqrt(ν0)/2π
    println(prefactor)
    plot!(l1_plot,β->exp(-ΔV₀₁*β)*prefactor,linestyle=styles[i],color=colors[i],label="",yaxis=:log)
end
savefig(l1_plot,"plots/l1.pdf")