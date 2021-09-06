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


using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm_functions_structs.jl")

################################################################################
#                                                      						   #
#  An individual-based model of host-parasite coevolution in continuous space  #
#                                                      						   #
#  Used to compare with results found under the SPDE model        			   #
#                                                      						   #
################################################################################

# parameter values
μₕ = 0.1;	# mutation rates
μₚ = 0.1;	
Eₕ = 0.01;	# environmental deviations
Eₚ = 0.01;
σₕ = 0.01;	# dispersal distances
σₚ = 0.01;
κₕ = 0.99;	# competition effects
κₚ = 0.99;
Rₕ = 0.02;	# competition radii
Rₚ = 0.02;
ιₕ = 0.95;	# interaction effects
ιₚ = 1.1;
Rᵢ = 0.02;	# interaction radius
πₘ = 0.9;	# max infection probability
γ = 0.05;	# infection sensitivity
αₕ = 1.4;	# abiotic effects
αₚ = 1.2;
Aₕ = 0.001;	# abiotic selection strengths
Aₚ = 0.001;
θ₀ₕ = 0;	# abiotic optima
θ₀ₚ = 0;

prs = hp_pars(μₕ=μₕ, μₚ=μₚ, Eₕ=Eₕ, Eₚ=Eₚ, σₕ=σₕ, σₚ=σₚ, θ₀ₕ=θ₀ₕ, θ₀ₚ=θ₀ₚ, κₕ=κₕ, κₚ=κₚ, 
	Rₕ=Rₕ, Rₚ=Rₚ, ιₕ=ιₕ, ιₚ=ιₚ, Rᵢ=Rᵢ, πₘ=πₘ, γ=γ, αₕ=αₕ, αₚ=αₚ, Aₕ=Aₕ, Aₚ=Aₚ);

# number of generations to halt at
T = 200;

# initial population sizes
nₕ = 100;
nₚ = 100;
n₀ = [nₕ nₚ];

# run the sim
Xₕ = sim(prs,n₀,T);

# pull out abundance time-series
nₕₕ = zeros(T);
nₚₕ = zeros(T);
for i in 1:T
	nₕₕ[i] = Xₕ[i].nₕ
	nₚₕ[i] = Xₕ[i].nₚ
end

# abundance dynamics
plot(1:T,nₚₕ)
plot(1:T,nₕₕ)

# spatial distribution in final generation
scatter(Xₕ[T].xₕ[:,1],Xₕ[T].xₕ[:,2])
scatter(Xₕ[T].xₚ[:,1],Xₕ[T].xₚ[:,2])

# trait distributions
histogram(Xₕ[T].zₕ)
histogram(Xₕ[T].zₚ)

# export space-trait data in final gen to a csv
host_df = DataFrame(id = 1:Xₕ[T].nₕ, trait = Xₕ[T].zₕ, x1 = Xₕ[T].xₕ[:,1], x2 = Xₕ[T].xₕ[:,2]);
para_df = DataFrame(id = 1:Xₕ[T].nₚ, trait = Xₕ[T].zₚ, x1 = Xₕ[T].xₚ[:,1], x2 = Xₕ[T].xₚ[:,2]);

CSV.write("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/host_df.csv",host_df)
CSV.write("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/para_df.csv",para_df)

# need to export multiple sets of results for use with ncf R pkg...