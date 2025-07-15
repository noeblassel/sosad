using DelimitedFiles,Statistics,ProgressMeter,Base.Threads

β = 1.677 # mol / kcal for T = 300 K

bins_sides_phi = -235:2:125.0
bin_centers_phi = -234:2:124
bin_sides_psi = -120:2:240
bin_centers_psi = -119:2:239

dϕ = bins_sides_phi.step
dψ = bin_centers_psi.step

binsize = Float64(dϕ*dψ)

get_bin(ϕ,ψ) = (1+floor(Int,(ϕ+235)/2),1+floor(Int,(ψ+120)/2))

datadir = "trajectories"
files = readdir(datadir)

nfiles = length(files)

observables = [:sum_a11_w,:sum_a22_w,:sum_a12_w,:sum_w,:var_a11_w,:var_a22_w,:var_a12_w,:var_w,:cov_a11_w_w,:cov_a22_w_w,:cov_a12_w_w,:nsamp]

hist = zeros(180,180,length(observables),nfiles)
burnin = 500
min_samp = 2 # to estimate variances and covariances

data_lock = ReentrantLock()

p = Progress(nfiles)

@threads for j=1:nfiles # one file per window
    f = files[j]
    M = Float64.(readdlm(joinpath(datadir,f),comments=true)[:,[1,2,3,5,7,9,11]])
    N = size(M,1)

    data = Dict{Tuple{Int,Int},Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}()

    for i=burnin:N
        binϕ,binψ = get_bin(M[i,2],M[i,3])
        if (binϕ,binψ)∈ keys(data)
            for (obs,vect)=zip(M[i,4:end],data[(binϕ,binψ)])
                push!(vect,obs)
            end
        else
            data[(binϕ,binψ)] = ([M[i,4]],[M[i,5]],[M[i,6]],[M[i,7]])
        end
    end

    max_vars = [-Inf,-Inf,-Inf,-Inf]

    for k = keys(data)
        a11,a22,a12,U = data[k]

        nsamp = length(U)

        w = exp.(β*U)
        a11_w = a11 .* w
        a22_w = a22 .* w
        a12_w = a12 .* w

        if length(U) >= min_samp
            var_a11_w = var(a11_w)
            var_a12_w = var(a12_w)
            var_a22_w = var(a22_w)

            var_w = var(w)

            cov_a11_w_w = cov(a11_w,w)
            cov_a12_w_w = cov(a12_w,w)
            cov_a22_w_w = cov(a22_w,w)

            max_vars = max.(max_vars,[var_a11_w,var_a22_w,var_a12_w,var_w])
            @lock data_lock hist[k...,5:11,j] .= [var_a11_w,var_a22_w,var_a12_w,var_w,cov_a11_w_w,cov_a22_w_w,cov_a12_w_w]
        else
            @lock data_lock hist[k...,5:8,j] .= 10000
            @lock data_lock hist[k...,9:11,j] .= 0
        end

        @lock data_lock hist[k...,1,j] = sum(a11_w)
        @lock data_lock hist[k...,2,j] = sum(a22_w)
        @lock data_lock hist[k...,3,j] = sum(a12_w)
        @lock data_lock hist[k...,4,j] = sum(w)

        @lock data_lock hist[k...,12,j] = nsamp
    end

    for l=1:4
        @lock data_lock replace!(hist[:,:,4+l,j],10000=>max_vars[l])
    end
    next!(p)
end

finish!(p)

a11_hist = zeros(180,180)
a22_hist = zeros(180,180)
a12_hist = zeros(180,180)

for i=1:180
    for j=1:180
        I = hist[i,j,12,:] .> 0 # indexes windows where samples are available

        a11_w = hist[i,j,1,I] ./ hist[i,j,12,I]
        a22_w = hist[i,j,2,I] ./ hist[i,j,12,I]
        a12_w = hist[i,j,3,I] ./ hist[i,j,12,I]
        w =  hist[i,j,4,I] ./ hist[i,j,12,I]

        a11 = a11_w ./ w
        a22 = a22_w ./ w
        a12 = a12_w ./ w

        var_a11_w = hist[i,j,5,I]
        var_a22_w = hist[i,j,6,I]
        var_a12_w = hist[i,j,7,I]
        var_w = hist[i,j,8,I]

        cov_a11_w_w = hist[i,j,9,I]
        cov_a22_w_w = hist[i,j,10,I]
        cov_a12_w_w = hist[i,j,11,I]

        #inversely weighted with delta method approximated variance -- **assumes constant correlation times across different windows at a given point --> poor results**

        # weight_a11 = @. inv((var_a11_w/w^2 + a11_w^2*var_w/w^4-2*a11_w*cov_a11_w_w/w^3) / hist[i,j,12,I])
        # weight_a12 = @. inv((var_a12_w/w^2 + a12_w^2*var_w/w^4-2*a12_w*cov_a12_w_w/w^3) / hist[i,j,12,I])
        # weight_a22 = @. inv((var_a22_w/w^2 + a22_w^2*var_w/w^4-2*a22_w*cov_a22_w_w/w^3) / hist[i,j,12,I])

        # inversely weighted by sample variance -- **assumes no correlation and constant correlation times --> better, but artifacts remain**
        # weight_a11 = inv.(hist[i,j,5,I] ./ hist[i,j,12,I])
        # weight_a22 = inv.(hist[i,j,6,I] ./ hist[i,j,12,I])
        # weight_a12 = inv.(hist[i,j,7,I] ./ hist[i,j,12,I])

        # weighted by number of samples -- **assumes same asymptotic variance in each window --> rather good results**
        # weight_a11 = hist[i,j,12,I]
        # weight_a22 = hist[i,j,12,I]
        # weight_a12 = hist[i,j,12,I]

        # uniform weighting --> simple and most effective
        weight_a11 = ones(size(a11))
        weight_a22 = ones(size(a22))
        weight_a12 = ones(size(a12))

        a11_hist[i,j] = sum(a11 .* weight_a11) / sum(weight_a11)
        a12_hist[i,j] = sum(a12 .* weight_a12) / sum(weight_a12)
        a22_hist[i,j] = sum(a22 .* weight_a22) / sum(weight_a22)
            
    end
end

writedlm("a11.out",a11_hist)
writedlm("a22.out",a22_hist)
writedlm("a12.out",a12_hist)
