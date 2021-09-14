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

include("ibm_functions_structs.jl")

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
σₕ = 0.001;	# dispersal distances
σₚ = 0.001;
κₕ = 0.95;	# competition effects
κₚ = 0.975;
Rₕ = 0.01;	# competition radii
Rₚ = 0.01;
ιₕ = 0.99;	# interaction effects
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
X1 = sim(prs,n₀,T);
X2 = sim(prs,n₀,T);
X3 = sim(prs,n₀,T);
X4 = sim(prs,n₀,T);
X5 = sim(prs,n₀,T);

# pull out abundance time-series
# nₕₕ = zeros(T);
# nₚₕ = zeros(T);
# for i in 1:T
# 	nₕₕ[i] = X1[i].nₕ
# 	nₚₕ[i] = X1[i].nₚ
# end

# abundance dynamics
# plot(1:T,nₚₕ)
# plot(1:T,nₕₕ)

# spatial distribution in final generation
# scatter(Xₕ[T].xₕ[:,1],Xₕ[T].xₕ[:,2])
# scatter(Xₕ[T].xₚ[:,1],Xₕ[T].xₚ[:,2])

# trait distributions
# histogram(Xₕ[T].zₕ)
# histogram(Xₕ[T].zₚ)

# export space-trait data in final gen to a csv
host_df1 = DataFrame(id = 1:X1[T].nₕ, trait = X1[T].zₕ, x1 = X1[T].xₕ[:,1], x2 = X1[T].xₕ[:,2]);
para_df1 = DataFrame(id = 1:X1[T].nₚ, trait = X1[T].zₚ, x1 = X1[T].xₚ[:,1], x2 = X1[T].xₚ[:,2]);
host_df2 = DataFrame(id = 1:X2[T].nₕ, trait = X2[T].zₕ, x1 = X2[T].xₕ[:,1], x2 = X2[T].xₕ[:,2]);
para_df2 = DataFrame(id = 1:X2[T].nₚ, trait = X2[T].zₚ, x1 = X2[T].xₚ[:,1], x2 = X2[T].xₚ[:,2]);
host_df3 = DataFrame(id = 1:X3[T].nₕ, trait = X3[T].zₕ, x1 = X3[T].xₕ[:,1], x2 = X3[T].xₕ[:,2]);
para_df3 = DataFrame(id = 1:X3[T].nₚ, trait = X3[T].zₚ, x1 = X3[T].xₚ[:,1], x2 = X3[T].xₚ[:,2]);
host_df4 = DataFrame(id = 1:X4[T].nₕ, trait = X4[T].zₕ, x1 = X4[T].xₕ[:,1], x2 = X4[T].xₕ[:,2]);
para_df4 = DataFrame(id = 1:X4[T].nₚ, trait = X4[T].zₚ, x1 = X4[T].xₚ[:,1], x2 = X4[T].xₚ[:,2]);
host_df5 = DataFrame(id = 1:X5[T].nₕ, trait = X5[T].zₕ, x1 = X5[T].xₕ[:,1], x2 = X5[T].xₕ[:,2]);
para_df5 = DataFrame(id = 1:X5[T].nₚ, trait = X5[T].zₚ, x1 = X5[T].xₚ[:,1], x2 = X5[T].xₚ[:,2]);

CSV.write("/ibm/host_df1.csv",host_df1)
CSV.write("/ibm/para_df1.csv",para_df1)
CSV.write("/ibm/host_df2.csv",host_df2)
CSV.write("/ibm/para_df2.csv",para_df2)
CSV.write("/ibm/host_df3.csv",host_df3)
CSV.write("/ibm/para_df3.csv",para_df3)
CSV.write("/ibm/host_df4.csv",host_df4)
CSV.write("/ibm/para_df4.csv",para_df4)
CSV.write("/ibm/host_df5.csv",host_df5)
CSV.write("/ibm/para_df5.csv",para_df5)

# need to export multiple sets of results for use with ncf R pkg...