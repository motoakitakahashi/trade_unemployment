# many country Krugman version of Helpman & Itskhoki (2010)
# Motoaki Takahashi
# May 2020


# packages
using LinearAlgebra, Statistics, Compat, Plots

# Parameters
N = 3 # number of countries
σ = 4.0 # elasticity of substitution
χ = 0.5 # labor force share in the Cobb-Douglas matching function
ζ = [1, 1, 1] # vacancy cost vector
f = [1, 1, 1] # entry cost vector

function VLratio(Φ, P, z)
    output = (σ - 1) ^ (1/χ) * (2*σ - 1) ^ (σ / (χ * (1-σ))) * Φ .^ (1 / ((σ-1)*χ) ) .* z .^ (1/χ) .* ζ .^ (-1/χ) .* f .^ (1/(χ * (1-σ))) .* P .^ (1/( χ * (1-σ) ))
    return output
end

function expenditure(θ, M, P)
    output = ζ .* θ .* L + M .* f .* P
    return output
end

function price_index(θ, M, z, t)
    output1 = (t') .^ (1 - σ) * (M .* ζ .^ (1 - σ) .* z .^ (σ - 1) .* θ .^ (χ * (1 - σ)))
    output2 = (2 * σ - 1) * (σ - 1) ^ (-1) * output1 .^ (1/(1-σ))
    return output2
end

function mass_of_firms(Φ, θ, z, L)
    output = (2*σ-1)^σ * (σ-1)^(-σ) * Φ .^ (-σ) .* ζ .^ σ .* θ .^ (χ * (σ-1) + 1) .* z .^ (1-σ) .* L
    return output
end

function market_potential(P, X, t)
    output1 = t .^ (1-σ) * (P .^ (σ - 1) .* X)
    output2 = output1 .^ (1/σ)
    return output2
end


# If we work in a while loop, we cannot access variabled defined outside the while loop.
# To deal with this issue, I wrap the while loop that computes an equilibrium with function.

# z is the productivity vector
# L is the labor force vector
function equilibrium(Φ_guess, P_guess, z, L, t)
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


    while dif > tol && it < maxit
        θ = VLratio(Φ, P, z)
        M = mass_of_firms(Φ, θ, z, L)
        X = expenditure(θ, M, P)
        P_new = price_index(θ, M, z, t)
        P_new = P_new ./ P_new[1] # the composite good in country 1 is the numeraire
        Φ_new = market_potential(P, X, t)


        dif = norm(([Φ_new; P_new] - [Φ; P]) ./ [Φ; P], Inf)
        Φ = λ * Φ_new + (1 - λ) * Φ
        P = λ * P_new + (1 - λ) * P
        it = it + 1
    end

    w = ζ .* θ .^ χ # nominal wage
    e = θ .^ (1-χ) # employment rate
    wel = e .* w ./ P # welfare (home production is zero)

    return [w./P; (ones(N, 1) - e); wel; dif]
end

#
# Comparative statistics w.r.t. productivity
z_grid = 1:0.01:3
result = zeros(3 * N + 1, length(z_grid))
L = [2, 2, 2]
t = [1 1.1 1.1; 1.1 1 1.1; 1.1 1.1 1] # trade cost matrix
# in t, the rows are origins and the columns are destinations

for i = 1:length(z_grid)
    z = [2, 2, z_grid[i]]
    result[:, i] = equilibrium([1, 1, 1], [1, 1, 1], z, L, t)
end
result

# plot comparative statics with respect to z

# real wages
plot(z_grid, result[1, :], label = "Country 1 and 2", legend = :topleft)
plot!(z_grid, result[3, :], label = "Country 3")
xlabel!("Productivity in Country 3")
ylabel!("Real Wage")

savefig("z_real_wage.pdf")
savefig("z_real_wage.png")

# unemployment
plot(z_grid, result[4, :], label = "Country 1 and 2", legend = :bottomleft)
plot!(z_grid, result[6, :], label = "Country 3")
xlabel!("Productivity in Country 3")
ylabel!("Unemployment Rate")

savefig("z_unemployment.pdf")
savefig("z_unemployment.png")

# welfare
plot(z_grid, result[7, :], label = "Country 1 and 2", legend = :topleft)
plot!(z_grid, result[9, :], label = "Country 3")
xlabel!("Productivity in Country 3")
ylabel!("Welfare")

savefig("z_welfare.pdf")
savefig("z_welfare.png")

#
# Comparative statics with respect to L
L_grid = 1 : 0.01 : 3
z = [2, 2, 2]
t = [1 1.1 1.1; 1.1 1 1.1; 1.1 1.1 1]
result = zeros(3 * N + 1, length(L_grid))

for i = 1:length(L_grid)
    L = [2, 2, L_grid[i]]
    result[:, i] = equilibrium([1, 1, 1], [1, 1, 1], z, L, t)
end
result

# real wages
plot(L_grid, result[1, :], label = "Country 1 and 2", legend = :topleft)
plot!(L_grid, result[3, :], label = "Country 3")
xlabel!("Labor in Country 3")
ylabel!("Real Wage")

savefig("L_real_wage.pdf")
savefig("L_real_wage.png")

# unemployment
plot(L_grid, result[4, :], label = "Country 1 and 2", legend = :topright)
plot!(L_grid, result[6, :], label = "Country 3")
xlabel!("Labor in Country 3")
ylabel!("Unemployment Rate")

savefig("L_unemployment.pdf")
savefig("L_unemployment.png")

# welfare
plot(L_grid, result[7, :], label = "Country 1 and 2", legend = :topleft)
plot!(L_grid, result[9, :], label = "Country 3")
xlabel!("Labor in Country 3")
ylabel!("Welfare")

savefig("L_welfare.pdf")
savefig("L_welfare.png")


# Comparatie statics with resoect to international trade costs
z = [2, 2, 2]
L = [2, 2, 2]

t_grid = 1:0.01:3
result = zeros(3 * N + 1, length(t_grid))

for i = 1:length(t_grid)
    inter_t = t_grid[i]
    t = repeat([inter_t], N, N) - Diagonal(fill(inter_t, N)) + I(N)
    result[:, i] = equilibrium([1, 1, 1], [1, 1, 1], z, L, t)
end
result

# plot real wages
plot(t_grid, result[1, :], label = "")
xlabel!("International Trade Cost")
ylabel!("Real Wage")

savefig("t_real_wage.pdf")
savefig("t_real_wage.png")

# plot unemployment
plot(t_grid, result[4, :], label = "")
xlabel!("International Trade Cost")
ylabel!("Unemployment Rate")

savefig("t_unemployment.pdf")
savefig("t_unemployment.png")

# plot welfare
plot(t_grid, result[7, :], label = "")
xlabel!("International Trade Cost")
ylabel!("Welfare")

savefig("t_welfare.pdf")
savefig("t_welfare.png")
