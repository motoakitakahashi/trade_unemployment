# many country Krugman version of Helpman & Itskhoki (2010)
# Motoaki Takahashi
# May 2020

# packages
using LinearAlgebra, Statistics, Compat

# Parameters
N = 3 # number of countries
σ = 4.0 # elasticity of substitution
χ = 0.5 # labor force share in the Cobb-Douglas matching function
z = [5, 5, 5] # productivity vector
ζ = [1, 1, 1] # vacancy cost vector
f = [1, 1, 1] # entry cost vector
t = [1 1.1 1.1; 1.1 1 1.1; 1.1 1.1 1] # trade cost matrix
# in t, the rows are origins and the columns are destinations
L = [1, 1, 1] # labor force vector

function VLratio(Φ, P)
    output = (σ - 1) ^ (1/χ) * (2*σ - 1) ^ (σ / (χ * (1-σ))) * Φ .^ (1/(σ-1)*χ) .* z .^ (1/χ) .* ζ .^ (-1/χ) .* f .^ (1/(χ * (1-σ))) .* P .^ (1/(χ * (1-σ)))
    return output
end

function expenditure(θ, M, P)
    output = ζ .* θ .* L + M .* f .* P
    return output
end

function price_index(θ, M)
    output1 = (2 * σ - 1) * (σ - 1) ^ (-1) * t' * (M .* ζ .* z .^ (-1) .* θ .^ χ)
    output2 = output1 .^ (1/(1-σ))
    return output2
end

function mass_of_firms(Φ, θ)
    output = (2*σ-1)^σ * (σ-1)^(-σ) * Φ .^ (-σ) .* ζ .^ σ .* θ .^ (χ * (σ-1) + 1) .* z .^ (1-σ) .* L
    return output
end

function market_potential(P, X)
    output1 = t .^ (1-σ) * (P .^ (σ - 1) .* X)
    output2 = output1 .^ (1/σ)
    return output2
end


# If we work in a while loop, we cannot access variables defined outside the while loop.
# To deal with this issue, I wrap the while loop that computes an equilibrium with function.

function equilibrium(Φ_guess, P_guess)
    # setup for the while loop
    maxit = 5000
    tol = 10 ^ (-6)
    it = 0
    λ = 0.5
    dif = 1

    Φ = Φ_guess
    P = P_guess
    θ = zeros(N, 1) # I need to write an object before a while loop...
    X = zeros(N, 1)
    M = zeros(N, 1)

    while dif > tol && it < maxit
        θ = VLratio(Φ, P)
        M = mass_of_firms(Φ, θ)
        X = expenditure(θ, M, P)
        P_new = price_index(θ, M)
        P_new = P_new ./ P_new[1] # the composite good in country 1 is the numeraire
        Φ_new = market_potential(P_new, X)


        dif = norm(([Φ_new; P_new] - [Φ; P]) ./ [Φ; P], Inf)
        Φ = λ * Φ_new + (1 - λ) * Φ
        P = λ * P_new + (1 - λ) * P
        it = it + 1
    end

    dif = fill(dif, N)
    w = ζ .* θ .^ χ
    u = ones(N, 1) - θ .^ (1-χ)

    return [Φ θ X P w u dif]
end

equilibrium([1, 1, 2], [1, 1, 1])
