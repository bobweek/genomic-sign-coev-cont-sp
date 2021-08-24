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

########################################################
#                                                      #
#  An individual-based model for the entire community  #
#                                                      #
#  Used to compare with diffusion approximation        #
#                                                      #
########################################################

# parameter values
μₕ = 0.1
μₚ = 0.1
Eₕ = 0.01
Eₚ = 0.01
σₕ = 0.01
σₚ = 0.01
κₕ = 0.99
κₚ = 0.99
Rₕ = 0.02
Rₚ = 0.02
ιₕ = 0.9
ιₚ = 1.1
Rᵢ = 0.02
πₘ = 0.9
γ = 0.1
αₕ = 1.2
αₚ = 1.1
Aₕ = 0.1
Aₚ = 0.1
θ₀ₕ =  1
θ₀ₚ = -1

# initial population sizes
nₕ = 1000
nₚ = 1000

# uniform random positions on unit square
xₕ = rand(nₕ,2)
xₚ = rand(nₚ,2)

# initial breeding values
gₕ = rand( Normal(θ₀ₕ,√μₕ), nₕ)
gₚ = rand( Normal(θ₀ₚ,√μₚ), nₚ)

# initial trait values
Eₕₘ = √Eₕ*Matrix(I, nₕ, nₕ)
zₕ = vec(rand(MvNormal(gₕ,Eₕₘ)))
Eₚₘ = √Eₚ*Matrix(I, nₚ, nₚ)
zₚ = vec(rand(MvNormal(gₚ,Eₚₘ)))

# set up initial population
X = hp_struct(zₕ=zₕ, zₚ=zₚ, gₕ=gₕ, gₚ=gₚ, nₕ=nₕ, nₚ=nₚ, xₕ=xₕ, xₚ=xₚ, μₕ=μₕ, μₚ=μₚ, Eₕ=Eₕ, Eₚ=Eₚ, σₕ=σₕ, σₚ=σₚ, 
	θ₀ₕ=θ₀ₕ, θ₀ₚ=θ₀ₚ, κₕ=κₕ, κₚ=κₚ, Rₕ=Rₕ, Rₚ=Rₚ, ιₕ=ιₕ, ιₚ=ιₚ, Rᵢ=Rᵢ, πₘ=πₘ, γ=γ, αₕ=αₕ, αₚ=αₚ, Aₕ=Aₕ, Aₚ=Aₚ)

# always a good idea to inspect a single iteration
X = update(X)

# number of generations to halt at
T = 10

# set up history of population
Xₕ = fill(X,T)

# simulate
for i in 2:T
	if Xₕ[i-1].nₕ>0 && Xₕ[i-1].nₚ>0
		Xₕ[i] = update(Xₕ[i-1])
	else
		Xₕ[i] = Xₕ[i-1]
	end
end


scatter(Xₕ[10].xₕ[:,1],Xₕ[10].xₕ[:,2])

scatter(Xₕ[10].xₚ[:,1],Xₕ[10].xₚ[:,2])
