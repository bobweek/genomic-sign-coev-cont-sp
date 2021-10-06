#
# script for checking how well theoretical prediction of model parameters holds up
#

# be sure to cd into the folder containing this file first
include("../ibm_functions_structs.jl")

# number of simulations to run per param combo
N = 3

# parameter values
μₕ = 0.1;# mutation rates
μₚ = 0.1;
Eₕ = 0.00;# environmental deviations
Eₚ = 0.00;
σₕ = 0.1;# dispersal distances
σₚ = 0.1;
κₕ = [0.75 0.8 0.975];# competition effects
κₚ = [0.85 0.9 0.985];
Rₕ = 0.001;# competition radii
Rₚ = 0.001;
ιₕ = 0.99;# interaction effects
ιₚ = 1.01;
Rᵢ = 0.01;# interaction radius
πₘ = 1.0;# max infection probability
γ = 0.1;# infection sensitivity
αₕ = 1.2;# abiotic effects
αₚ = 1.2;
Aₕ = 0.01;# abiotic selection strengths
Aₚ = 0.01;
θ₀ₕ = 0;# abiotic optima
θ₀ₚ = 0;

# values of competition effects to check
kₕ = 0.985
kₚ = 0.985

# solving for expected densities
function ρsolve(X)
    ρ̃ₕ = exp(X[1])
    ρ̃ₚ = exp(X[2])

    B̃ₕ = Bₕ(γ, πₘ, ιₚ, ρ̃ₚ, ρ̃ₕ)
    B̃ₚ = Bₚ(γ, πₘ, ιₕ, Rᵢ, ρ̃ₕ)
    ṽₚ = vₚ(Eₚ, μₚ, Aₚ, B̃ₚ)
    r̃ₕ = rₕ(αₕ, πₘ, ιₕ, γ, ṽₚ, ρ̃ₚ, ρ̃ₕ)
    r̃ₚ = rₚ(αₚ, πₘ, ιₚ, Rᵢ, ρ̃ₕ)

    val = √((ρ̃ₕ - ρₕ(kₕ, r̃ₕ, μₕ, Aₕ, B̃ₕ))^2 + (ρ̃ₚ - ρₚ(kₚ, r̃ₚ, μₚ, Aₚ, B̃ₚ))^2)

    return val

end

# # for checking exp densities, coll vars and char lngths
initp = rand(Uniform(1, 20), 2)
res = optimize(ρsolve, initp, SimulatedAnnealing(), Optim.Options(iterations = Int64(1e5)))
ρ̃ₕ = exp(Optim.minimizer(res)[1])
ρ̃ₚ = exp(Optim.minimizer(res)[2])
B̃ₕ = Bₕ(γ, πₘ, ιₕ, ρ̃ₚ, ρ̃ₕ)
B̃ₚ = Bₚ(γ, πₘ, ιₚ, Rᵢ, ρ̃ₕ)
ṽₚ = vₚ(Eₚ, μₚ, Aₚ, B̃ₚ)
ṽₕ = vₕ(Eₕ, μₕ, Aₕ, B̃ₕ)
r̃ₕ = rₕ(αₕ, πₘ, ιₕ, γ, ṽₚ, ρ̃ₚ, ρ̃ₕ)
r̃ₚ = rₚ(αₚ, πₘ, ιₚ, Rᵢ, ρ̃ₕ)
Vₕ(ρ̃ₕ, σₕ, Aₕ, B̃ₕ)
Vₚ(ρ̃ₚ, σₚ, Aₚ, B̃ₚ)
ξₕ(ṽₕ, σₕ, Aₕ, B̃ₕ)
ξₚ(ṽₚ, σₚ, Aₚ, B̃ₚ)

ncombo = length(κₕ) * length(κₚ)

prs = hp_pars(
    μₕ = μₕ,
    μₚ = μₚ,
    Eₕ = Eₕ,
    Eₚ = Eₚ,
    σₕ = σₕ,
    σₚ = σₚ,
    θ₀ₕ = θ₀ₕ,
    θ₀ₚ = θ₀ₚ,
    κₕ = κₕ[1],
    κₚ = κₚ[1],
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

nₕ = 100
nₚ = 100
n₀ = [nₕ nₚ]

# create an array to hold output from each simulation
X = fill(sim(prs, n₀, 1), ncombo)

cnt = 1
for i = 1:length(κₕ)
    for j = 1:length(κₚ)

        parm_df = DataFrame(
            μₕ = μₕ,
            μₚ = μₚ,
            Eₕ = Eₕ,
            Eₚ = Eₚ,
            σₕ = σₕ,
            σₚ = σₚ,
            θ₀ₕ = θ₀ₕ,
            θ₀ₚ = θ₀ₚ,
            κₕ = κₕ[i],
            κₚ = κₚ[j],
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
        CSV.write(string("parm_df", cnt, "-", k, ".csv"), parm_df)

        for k = 1:N

            prs = hp_pars(
                μₕ = μₕ,
                μₚ = μₚ,
                Eₕ = Eₕ,
                Eₚ = Eₚ,
                σₕ = σₕ,
                σₚ = σₚ,
                θ₀ₕ = θ₀ₕ,
                θ₀ₚ = θ₀ₚ,
                κₕ = κₕ[i],
                κₚ = κₚ[j],
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
            nₕ = 100
            nₚ = 100
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

            CSV.write(string("host_df", cnt, "-", k, ".csv"), host_df)
            CSV.write(string("para_df", cnt, "-", k, ".csv"), para_df)

            cnt += 1
        end
    end
end


