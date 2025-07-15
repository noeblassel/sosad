k = 0.01293 # linear factor
s = 4 # scale factor
e = 0.7 # energy scale

V(x) = e*(1-cos(s*x) -exp(-(s*x-1)^2/2)+k*s*x)
V′(x) = e*(s*(sin(s*x) + (s*x-1)*exp(-(s*x-1)^2/2)))
V′′(x) = e*(s^2*(cos(s*x)-s*x*(s*x-2)*exp(-(s*x-1)^2/2)))

mu(x,β) = exp(-β*V(x))

# plot(V,0,1,label="V")

N = 100000
X = collect(range(-6/s,6/s,N+1))
pop!(X)

@inline prev(i) = mod(i-1,1:N)
@inline next(i) = mod(i+1,1:N)


minima = [i for i=2:N-1 if  V(X[prev(i)])>V(X[i])<V(X[next(i)])]
saddles = [i for i=2:N-1 if V(X[prev(i)])<V(X[i])>V(X[next(i)])]

i0 = first(minima)
z0,z1 = saddles

ΔV₀₁ = min(V(X[z0])-V(X[i0]),V(X[z1])-V(X[i0]))
# ΔV₁₀ = min(V(X[z0])-V(X[i1]),V(X[z1])-V(X[i1]))

# lowest_saddle = (V(X[z0]) < V(X[z1])) ? z0 : z1
# highest_saddle = (lowest_saddle == z0) ? z1 : z0

ν0 = V′′(X[i0])
ν1 = -V′′(X[z0])
ν2 = -V′′(X[z1])

# ℓ = min(κ0,κ1,V′′(X[i0]))
# κ_high = -V′′(X[highest_saddle])
# κ_low = -V′′(X[lowest_saddle])

# V_high = V(X[highest_saddle])
# V_low = V(X[lowest_saddle])

#Ωind = z0:z1 # bassin of i0
#Ωind = mod.(z1-N:z0,(1:N,)) # bassin of i1
# Ωind = [i for i=1:N if (X[z1]+0.1>X[i])&&(X[i]>X[z0]-0.1)] # domain encompassing bassin of i0
# Ω = X[Ωind]


println("---- Critical points of V ---")
println("\t--- Minima")
println("\t\t x0: ", X[i0],"\n")
println("\t--- Saddles")
println("\t\t z0: ", X[z0])
println("\t\t z1: ", X[z1],"\n")

println("--- Values of V ---")

println("\t--- Minima")
println("\t\t V(x0): ", V(X[i0]),"\n")
# println("\t\t x1: ", V(X[i1]),"\n")

println("\t--- Saddles")
println("\t\t V(z0): ", V(X[z0]))
println("\t\t V(z1): ", V(X[z1]),"\n")

println("\t--- Potential barriers ---")
println("\t\t x0 → x1: ", ΔV₀₁,"\n")
# println("\t\t x1 → x0: ", ΔV₁₀,"\n")


println("--- Values of V′′")

println("\t--- Minima")
println("\t\t V′′(x0): ", V′′(X[i0]),"\n")
# println("\t\t x1: ", V′′(X[i1]),"\n")

println("\t--- Saddles")
println("\t\t V′′(z0): ", V′′(X[z0]))
println("\t\t V′′(z1): ", V′′(X[z1]))

Ω_f(β,a,b) = [i for i=1:N if X[z0]-a/sqrt(β)<X[i]<X[z1]+b/sqrt(β)]