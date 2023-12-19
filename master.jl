using ProgressMeter, Plots, StatsPlots, StatsBase, BivariateCopulas, Dates, DataFrames
using SQLite, DBInterface
using Turing, BivariateCopulas, EmpiricalCopulas, KernelDensity

#Basis for the test statistics
function Qₙ(x; y)
    x = typeof(x) <: Real ? [x] : x
    κ = length(x)
    Δy = [0.,  diff(y)...]
 
    n = length(y)
    Q = 0.
    
    for t in (κ+1):n
        Q += Δy[t] * prod(y[t-k] ≤ x[k] for k in 1:κ)
    end
 
    1/sqrt(n) * Q
end

#Cramer-von Mises
function Tₙ(y)
    n = length(y)
 
    1/(n-1) * sum(Qₙ(y[t-1]; y = y)^2 for t in 2:n)
end

#Kolmogorov-Smirnov
function Sₙ(y)
    S = map(z -> Qₙ(z; y = y), y)
    return maximum(abs.(S))
end
 
𝐒 = Float64[]
𝐓 = Float64[]
 
n = 7500

@showprogress for _ in 1:n
    y = cumsum(randn(n))

    push!(𝐓, Tₙ(y))
    push!(𝐒, Sₙ(y))
end

#Density of the two test statistics
StatsPlots.density(𝐒, label = "Kolmogorov-Smirnov")
StatsPlots.density!(𝐓, label = "Cramer-von Mises")

#Critical values of Kolmogorov-Smirnov and Cramer-von Mises
S₀ = quantile(𝐓, [0.01, 0.05, 0.1, 0.9, 0.95, 0.99])
T₀ = quantile(𝐒, [0.01, 0.05, 0.1, 0.9, 0.95, 0.99])

#Testing on data
Results = []

@showprogress for _ in _
    sort!(_, :)

    x = _[:, :]

    Δx = diff(x)

    m = mean(Δx)
    s = std(Δx .- m)
    n = length(Δx)

    #Testing if increments have mean zero 
    ℋ₀ = m - 1.96*(s/sqrt(n)) ≤ 0 ≤ m + 1.96*(s/sqrt(n))   

    #Standize 
    Δz = Δx / s
    z = cumsum(Δz)

    #Find test statistics
    S = length(z) > 0 ? Sₙ(z) : Inf
    T = length(z) > 0 ? Tₙ(z) : Inf 

    ℋ₁ = 0 .≤ S .≤ S₀
    ℋ₂ = 0 .≤ T .≤ T₀

    push!(Results, (m, s, ℋ₀, ℋ₁, ℋ₂, n, S, T))
end

#Plot simulated against real data for Kolmogorov-Smirnov and Cramer-von Mises
𝐒′ = [r[7] for r in Results]
StatsPlots.histogram(𝐒, label = "Simulation", normalize = true)
StatsPlots.histogram!(𝐒′, label = "Data", normalize = true)
StatsPlots.density!(𝐒, label = "Simulation (PDF)")
StatsPlots.density!([s for s in 𝐒′ if 0 ≤ t < Inf], label = "Data (PDF)")

𝐓′ = [r[8] for r in Results]
StatsPlots.histogram(𝐓, label = "Simulation", normalize = true)
StatsPlots.histogram!(𝐓′, label = "Data", normalize = true)
StatsPlots.density!(𝐓, label = "Simulation (PDF)")
StatsPlots.density!([t for t in 𝐓′ if 0 ≤ t < Inf], label = "Data (PDF)")

#Critical values of Kolmogorov-Smirnov and Cramer-von Mises and if increments have mean zero
sum([r[3] for r in Results])/length(Results) #Increments mean zero
sum([r[4] for r in Results])/length(Results) #Kolmogorov-Smirnov
sum([r[5] for r in Results])/length(Results) #Cramer-von Mises

#Testing if the difference has mean zero
Results2 = []
for _ in _
    Δ = _[:, :]
    n = length(Δ)
    m = mean(Δ)
    s = std(Δ)
    ℋ₃ = m - 1.96*(s/sqrt(n)) ≤ 0 ≤ m + 1.96*(s/sqrt(n))

    push!(Results2, ℋ₃)
end

#Check how many times the test is rejected
StatsPlots.scatter(Results2) 
sum(Results2)/length(Results2) 

#Copula fitting
time_to_delivery = sort(unique([:, :]))
l = 8
data = filter(: => x -> x == time_to_delivery[l], :)

X = data[:, :]
Y = data[:, :]

Fₙ = ecdf(X)
Gₙ = ecdf(Y)
fₙ(x) = pdf(kde(X),x)
gₙ(y) = pdf(kde(Y),y)

Û = Fₙ.(X)
V̂ = Gₙ.(Y)

ĉ = BernsteinCopula(Û,V̂,10) 

z = [pdf(ĉ, [Fₙ(x),Gₙ(y)])*fₙ(x)*gₙ(y) for x in LinRange(extrema(X)...,100), y in LinRange(extrema(Y)...,100)]
heatmap(LinRange(extrema(X)...,100), LinRange(extrema(Y)...,100), z')
scatter!(X,Y, alpha = 0.15, label = "Data points")

n = length(data[:, 1])
k = 50
D = 1:k:n

indices = Int64[]

for i in 1:n 
    push!(indices, searchsortedlast(D,i))
end

data[:, :] = indices

fₙ(x) = pdf(kde(X),x)
K(x,y) = pdf(ĉ, [Fₙ.(x),Gₙ.(y)])*fₙ(x)

i = LinRange(extrema(X)...,100)
pieces = LinRange(minimum(Y),maximum(Y),3)

plot()
for y in pieces
    plot!(i, x -> K(x,y), label="y = $y")
end
title!("")
