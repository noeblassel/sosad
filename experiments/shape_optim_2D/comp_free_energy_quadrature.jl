using Plots,Cubature,DelimitedFiles,ProgressMeter,Zygote


function main(ε,β,nbins)
    V(x,y) = (x^2-1)^2 + (x^2+y^2-1)^2/ε + 1/sqrt(x^2+y^2)
    V(q) = V(q...)

    ξ1(x,y) = atan(y,x)/π
    ξ1(q) = ξ1(q...)
    sq∇ξ1(q) = sum(abs2,first(gradient(ξ1,q)))

    # ξ2(x,y) = x
    # sq∇ξ2(x,y) = 1.0

    ξ1min, ξ1max = -1.0,1.0
    ξ1_perp_min,ξ1_perp_max = 0.0, 5.0

    ξ2min,ξ2max= -2.0,2.0
    ξ2_perp_min,ξ2_perp_max = -5.0, 5.0

    F1 = zeros(nbins)
    D1 = zeros(nbins)

    F2 = zeros(nbins)
    D2 = ones(nbins) # analytical

    h1 = (ξ1max -ξ1min)/nbins
    h2 = (ξ2max - ξ2min)/nbins

    ξ1s = collect(range(ξ1min+h1/2,ξ1max-h1/2,nbins))
    ξ2s = collect(range(ξ2min+h2/2,ξ2max-h2/2,nbins))

    max_err_F1 = -Inf
    max_err_D1 = -Inf
    max_err_F2 = -Inf

    @showprogress for k=1:nbins
        σ1(t) = [t*cos(π*ξ1s[k]),t*sin(π*ξ1s[k])] # unit parametrization of ξ1 level set
        σ2(t) = [ξ2s[k],t] # unit parametrization of ξ2 level set

        Z1,err = hquadrature(t-> begin q = σ1(t); return exp(-β * V(q))/sqrt(sq∇ξ1(q)) end,ξ1_perp_min,ξ1_perp_max) # local partition function
        F1[k] = -log(Z1)/β # free energy
        max_err_F1 = max(max_err_F1,err)
        D1[k],err = hquadrature(t-> begin q = σ1(t); return exp(-β * V(q))*sqrt(sq∇ξ1(q))/Z1 end,ξ1_perp_min,ξ1_perp_max) # effective diffusion coeff
        max_err_D1 = max(max_err_D1,err)

        Z2,err = hquadrature(t-> begin q = σ2(t); return exp(-β * V(q)) end,ξ2_perp_min,ξ2_perp_max)
        F2[k] = -log(Z2)/β
        max_err_F2 = max(max_err_F2,err)
    end

    M = [ξ1s (F1 .- minimum(F1)) D1 ξ2s (F2 .- minimum(F2)) D2]

    writedlm("./series/FD_quad_eps=$(ε)_beta=$(β).out",M)

    return max_err_F1,max_err_D1,max_err_F2
end


εs = [2.0,1.0,0.5,0.1]
β = 2.0

for ε = εs
    println(ε)
    main(ε,2.0,1000)
end
