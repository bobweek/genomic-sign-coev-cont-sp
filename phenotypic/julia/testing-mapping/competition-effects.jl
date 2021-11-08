#
# script for checking how well diffusion-limit prediction of 
# density and additive genetic variance holds up to simulated data
#

# be sure to cd into the folder containing this file first
include("../ibm_functions_structs.jl")

# number of simulations to run per parameter combo
n = 3

# parameter values
μ = [0.1 0.5 1.0]; # mutation rates
Eₕ = 0.00;  # environmental deviations
Eₚ = 0.00;
σₕ = 0.01;   # dispersal distances
σₚ = 0.01;
κ = [0.75 0.8 0.9]; # competition effects
Rₕ = 0.01; # competition radii
Rₚ = 0.01;
ιₕ = 0.99;  # interaction effects
ιₚ = 1.01;
Rᵢ = 0.01;  # interaction radius
πₘ = 1.0;   # max infection probability
γ = 0.1;    # infection sensitivity
αₕ = 1.175;   # abiotic effects
αₚ = 1.1575;
Aₕ = 0.02;  # abiotic selection strengths
Aₚ = 0.02;
θ₀ₕ = 0.0;  # abiotic optima
θ₀ₚ = 0.0;

ncombo = length(κ) * length(μ)

prs = hp_pars(
    μₕ = μ[1],
    μₚ = μ[1],
    Eₕ = Eₕ,
    Eₚ = Eₚ,
    σₕ = σₕ,
    σₚ = σₚ,
    θ₀ₕ = θ₀ₕ,
    θ₀ₚ = θ₀ₚ,
    κₕ = κ[3],
    κₚ = κ[3],
    Rₕ = Rₕ,
    Rₚ = Rₚ,
    ιₕ = ιₕ,
    ιₚ = ιₚ,
    Rᵢ = Rᵢ,
    πₘ = πₘ,
    γ = γ,
    αₕ = αₕ,
    αₚ = αₚ,
    Aₕ = Aₕ,
    Aₚ = Aₚ,
)

nₕ = 1000
nₚ = 1000
n₀ = [nₕ nₚ]

# try a single simulation
T = 1000
Y = sim(prs, n₀, T, true)

histnₚ = zeros(T)
histnₕ = zeros(T)
for i in 1:T
    histnₚ[i] = Y[i].nₚ
    histnₕ[i] = Y[i].nₕ
end

plot(1:T,histnₚ)
plot(1:T,histnₕ)

# create an array to hold output from each simulation
X = fill(sim(prs, n₀, 1), ncombo)

cnt = 1
for i = 1:length(κ)
    for j = 1:length(μ)

        parm_df = DataFrame(
            μₕ = μ[j],
            μₚ = μ[j],
            Eₕ = Eₕ,
            Eₚ = Eₚ,
            σₕ = σₕ,
            σₚ = σₚ,
            θ₀ₕ = θ₀ₕ,
            θ₀ₚ = θ₀ₚ,
            κₕ = κ[i],
            κₚ = κ[i],
            Rₕ = Rₕ,
            Rₚ = Rₚ,
            ιₕ = ιₕ,
            ιₚ = ιₚ,
            Rᵢ = Rᵢ,
            πₘ = πₘ,
            γ = γ,
            αₕ = αₕ,
            αₚ = αₚ,
            Aₕ = Aₕ,
            Aₚ = Aₚ,
        )
        CSV.write(string("parm_df-", cnt, ".csv"), parm_df)

        for m = 1:n

            prs = hp_pars(
                μₕ = μ[j],
                μₚ = μ[j],
                Eₕ = Eₕ,
                Eₚ = Eₚ,
                σₕ = σₕ,
                σₚ = σₚ,
                θ₀ₕ = θ₀ₕ,
                θ₀ₚ = θ₀ₚ,
                κₕ = κ[i],
                κₚ = κ[i],
                Rₕ = Rₕ,
                Rₚ = Rₚ,
                ιₕ = ιₕ,
                ιₚ = ιₚ,
                Rᵢ = Rᵢ,
                πₘ = πₘ,
                γ = γ,
                αₕ = αₕ,
                αₚ = αₚ,
                Aₕ = Aₕ,
                Aₚ = Aₚ,
            )

            # number of generations to halt at
            T = 1000

            # initial population sizes
            nₕ = 1000
            nₚ = 1000
            n₀ = [nₕ nₚ]

            # run the sim
            X[cnt] = sim(prs, n₀, T)

            host_df = DataFrame(
                id = 1:X[cnt][T].nₕ,
                trait = X[cnt][T].zₕ,
                x1 = X[cnt][T].xₕ[:, 1],
                x2 = X[cnt][T].xₕ[:, 2],
            )
            para_df = DataFrame(
                id = 1:X[cnt][T].nₚ,
                trait = X[cnt][T].zₚ,
                x1 = X[cnt][T].xₚ[:, 1],
                x2 = X[cnt][T].xₚ[:, 2],
            )

            CSV.write(string("host_df-", cnt, "-", m, ".csv"), host_df)
            CSV.write(string("para_df-", cnt, "-", m, ".csv"), para_df)

        end

        cnt += 1

    end
end
