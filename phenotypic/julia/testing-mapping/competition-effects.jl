################################################################################
##
## AUTHOR: Bob Week
##
## DATE: 08/23/2021
##
## In this script we simulate an individual-based model for hosts and parasites
## coevolving in continous space.
##
## This script depends on another script called "ibm_functions_structs.jl".
## In that script we provide definitions of data structures for model parameters
## and for state variables. In that script we also define methods for iterating
## the simulation. With all the gory details located elsewhere, we can focus
## this script on simulating the model for a specified duration.
##
################################################################################


using Parameters,
    Statistics,
    Random,
    LinearAlgebra,
    Distributions,
    StatsBase,
    StatsPlots,
    Plots,
    DataFrames,
    CSV,
    Optim

include("ibm_functions_structs.jl")

################################################################################
#                                                      						   #
#  An individual-based model of host-parasite coevolution in continuous space  #
#                                                      						   #
#  Used to compare with results found under the SPDE model        			   #
#                                                      						   #
################################################################################

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

for i = 1:3
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
        κₚ = κₚ[i],
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


