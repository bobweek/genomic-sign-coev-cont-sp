#
# script for checking how well theoretical prediction of model parameters holds up
#

include("../ibm_functions_structs.jl")

# number of simulations to run per param combo
N = 10

# parameter values
μₕ = 0.1;# mutation rates
μₚ = 0.1;
Eₕ = 0.001;# environmental deviations
Eₚ = 0.001;
σₕ = 0.001;# dispersal distances
σₚ = 0.001;
κₕ = [0.975 0.985 0.995];# competition effects
κₚ = [0.98 0.99 0.999];
Rₕ = 0.005;# competition radii
Rₚ = 0.005;
ιₕ = 0.99;# interaction effects
ιₚ = 1.1;
Rᵢ = 0.01;# interaction radius
πₘ = 1.0;# max infection probability
γ = 0.1;# infection sensitivity
αₕ = 1.4;# abiotic effects
αₚ = 1.2;
Aₕ = 0.002;# abiotic selection strengths
Aₚ = 0.002;
θ₀ₕ = 0;# abiotic optima
θ₀ₚ = 0;

kₕ = κₕ[1]
kₚ = κₚ[1]

# solving for expected densities
function ρsolve(X)
    ρ̃ₕ = X[1]
    ρ̃ₚ = X[2]
    
    B̃ₕ = Bₕ(γ,πₘ,ιₚ,Rᵢ,ρ̃ₕ)
    B̃ₚ = Bₚ(γ,πₘ,ιₕ,ρ̃ₚ,ρ̃ₕ)
    ṽₚ = vₚ(Eₚ,μₚ,Aₚ,B̃ₚ)
    r̃ₕ = rₕ(αₕ,πₘ,ιₕ,γ,ṽₚ,ρ̃ₚ,ρ̃ₕ)
    r̃ₚ = rₚ(αₚ,πₘ,ιₚ,Rᵢ,ρ̃ₕ)

    val = √((ρ̃ₕ-ρₕ(kₕ,r̃ₕ,μₕ,Aₕ,B̃ₕ))^2 + (ρ̃ₚ-ρₚ(kₚ,r̃ₚ,μₚ,Aₚ,B̃ₚ))^2)

    return val

end

initp = rand(Uniform(1,20), 2)
res = optimize(ρsolve, initp)
Optim.minimizer(res)

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


X = fill(sim(prs, n₀, 1), 3)

for i = 1:length(κₕ)
    for j = 1:length(κₚ)
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
        X[i] = sim(prs, n₀, T)

        host_df = DataFrame(
            id = 1:X[i][T].nₕ,
            trait = X1[T].zₕ,
            x1 = X1[T].xₕ[:, 1],
            x2 = X1[T].xₕ[:, 2],
        )
        para_df = DataFrame(
            id = 1:X[i][T].nₚ,
            trait = X1[T].zₚ,
            x1 = X1[T].xₚ[:, 1],
            x2 = X1[T].xₚ[:, 2],
        )

        CSV.write(string("ibm/test-mapping/host_df", i, ".csv"), host_df)
        CSV.write(string("ibm/test-mapping/para_df", i, ".csv"), para_df)

    end
end


