# many country Krugman version of Helpman & Itskhoki (2010)
# Motoaki Takahashi
# May 2020

# packages
using LinearAlgebra, Statistics, Compat

# Parameters
N = 3 # number of countries
σ = 4.0 # elasticity of substitution
χ = 0.5 # labor force share in the Cobb-Douglas matching function
ζ = [1, 1, 1] # vacancy cost vector
f = [1, 1, 1] # entry cost vector
t = [1 1.1 1.1; 1.1 1 1.1; 1.1 1.1 1] # trade cost matrix
# in t, the rows are origins and the columns are destinations

function VLratio(Φ, z)
    output = (σ - 1) ^ (1/χ) * (2*σ - 1) ^ (σ / (χ * (1-σ))) * Φ .^ (1/(σ-1)*χ) .* z .^ (1/χ) .* ζ .^ (-1/χ) .* f .^ (1/(χ * (1-σ)))
    return output
end

function expenditure(θ, L)
    output = ζ .* θ .* L
    return output
end

function price_index(Φ, θ, z, L)
    output1 = (σ-1) .^ (-(σ+1)) * (2*σ-1) .^ (σ+1) * t' * (Φ .^ (-σ) .* ζ .^ (σ+1) .* z .^ (-σ) .* θ .^ (χ*σ+1) .* L)
    output2 = output1 .^ (1/(1-σ))
    return output2
end

function market_potential(P, X)
    output1 = t' * (P .^ (σ - 1) .* X)
    output2 = output1 .^ (1/σ)
    return output2
end


# If we work in a while loop, we cannot access variabled defined outside the while loop.
# To deal with this issue, I wrap the while loop that computes an equilibrium with function.

# z is the productivity vector
# L is the labor force vector
function equilibrium(Φ_guess, z, L)
    # setup for the while loop
    maxit = 5000
    tol = 10 ^ (-6)
    it = 0
    λ = 0.5
    dif = 1

    Φ = Φ_guess
    θ = zeros(N, 1) # I need to write an object before a while loop...
    X = zeros(N, 1)
    P = zeros(N, 1)

    while dif > tol && it < maxit
        θ = VLratio(Φ, z)
        X = expenditure(θ, L)
        P = price_index(Φ, θ, z, L)
        P = P ./ P[1] # the composite good in country 1 is the numeraire
        Φ_new = market_potential(P, X)

        dif = norm((Φ_new - Φ)./Φ, Inf)
        Φ = λ * Φ_new + (1 - λ) * Φ
        it = it + 1
    end

    w = ζ .* θ .^ χ # nominal wage
    e = θ .^ (1-χ) # employment rate
    wel = e .* w ./ P # welfare (home production is zero)

    return [w./P; (ones(N, 1) - e); wel; dif]
end

# Comparative statistics w.r.t. productivity
z_grid = 1:0.01:3
result = zeros(3 * N + 1, length(z_grid))
L = [1, 1, 1]

for i = 1:length(z_grid)
    z = [2, 2, z_grid[i]]
    result[:, i] = equilibrium([1, 1, 1], z, L)
end
result
