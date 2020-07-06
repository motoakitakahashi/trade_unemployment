# Motoaki Takahashi

# packages
using LinearAlgebra, Statistics, Compat, Plots, SpecialFunctions

# Parameters
N = 2 # number of locations
θ = 4 # trade elasticity
σ = 4 # demand elasticity of substitution
G = gamma(1 + (1 - σ)/θ) ^ (1/(1-σ)) # the constant in the price index, from EK
β = 0.5 # labor share
ν = [1, 1] # vacancy cost in terms of final goods

T = [1, 1] # absolute advantage
b = [0.5, 0.5] # replacement rate
d = [1 1.2;
    1.2 1] # trade costs
# in d, the rows are exporters and the columns are importers. Opposite to EK.

function Φ(w, P)
    output1 = T .^ (θ) .* (fill(2, N) - b) .^ (-θ*β) .* w .^ (-θ*β) .* P .^ (-θ*(1-β))
    output2 = d' .^ (-θ) * output1
    return output2
end


function price(w, P)
    output = G * Φ(w, P) .^ (-1/θ)
    return output
end

function trade_share(w, P)
    suboutput1 = T .^ θ .* (fill(2, N) - b) .^ (-θ*β) .* w .^ (-θ*β) .* P .^ (-θ*(1-β))
    numer = repeat(suboutput1, 1, N) .* d .^ (-θ)
    denom = repeat(Φ(w, P)', N, 1)
    output = numer ./ denom
end

function tax(w, E, U)
    numer = sum(b .* w .* U)
    denom = sum(w .* E + b .* w .* U)
    output = numer / denom
    return output
end



repeat([2, 4], 1, 4)
