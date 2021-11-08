################################################################################
##
## AUTHOR: Bob Week ⚧
##
## DATE: 08/23/2021
##
## In this script we provide definitions of data structures for model parameters
## and for state variables used in our simulations of individual-based models.
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
    Optim,
    KernelDensity, # i don't think this one's needed
    NearestNeighbors

# data type that holds population and parameters

@with_kw mutable struct hp_struct
    zₕ::Vector{Float64}  # host trait values
    zₚ::Vector{Float64}  # parasite trait values
    gₕ::Vector{Float64}  # host breeding values
    gₚ::Vector{Float64}  # parasite breeding values
    xₕ::Matrix{Float64}  # host locations
    xₚ::Matrix{Float64}  # parasite locations
    nₕ::Int64            # number of host individuals
    nₚ::Int64            # number of parasite individuals
    # z̄ₕ::Float64        # host mean traits on 100x100 grid
    # z̄ₚ::Float64        # parasite mean traits on 100x100 grid
    # vₕ::Float64        # host phenotypic variances on 100x100 grid
    # vₚ::Float64        # parasite phenotypic variances on 100x100 grid
    # Gₕ::Float64        # host additive genetic variances on 100x100 grid
    # Gₚ::Float64        # parasite additive genetic variances on 100x100 grid
    μₕ::Float64          # variance of mutation for host
    μₚ::Float64          # variance of mutation for parasite
    Eₕ::Float64          # variance of environmental deviation for host
    Eₚ::Float64          # variance of environmental deviation for parasite
    σₕ::Float64          # dispersal distance of host
    σₚ::Float64          # dispersal distance of parasite
    κₕ::Float64          # fitness effect of spatial competition for host
    κₚ::Float64          # fitness effect of spatial competition for parasite
    Rₕ::Float64          # radius of spatial competition for host
    Rₚ::Float64          # radius of spatial competition for parasite
    ιₕ::Float64          # fitness effect of interspp interaction for host
    ιₚ::Float64          # fitness effect of interspp interaction for parasite
    Rᵢ::Float64          # radius of interspp interactions
    πₘ::Float64          # maximum probability of infection
    γ::Float64           # rate of decay of infection prb with trait difference
    αₕ::Float64          # maximum fitness effect of abiotic selection for host
    αₚ::Float64          # maximum fitness effect of abiotic selection for parasite
    Aₕ::Float64          # sensitivity of abiotic selection to trait value for host
    Aₚ::Float64          # sensitivity of abiotic selection to trait value for parasite
    θ₀ₕ::Float64         # baseline abiotic optimal trait value for host
    θ₀ₚ::Float64         # baseline abiotic optimal trait value for parasite
end

@with_kw mutable struct hp_pars
    μₕ::Float64          # variance of mutation for host
    μₚ::Float64          # variance of mutation for parasite
    Eₕ::Float64          # variance of environmental deviation for host
    Eₚ::Float64          # variance of environmental deviation for parasite
    σₕ::Float64          # dispersal distance of host
    σₚ::Float64          # dispersal distance of parasite
    κₕ::Float64          # fitness effect of spatial competition for host
    κₚ::Float64          # fitness effect of spatial competition for parasite
    Rₕ::Float64          # radius of spatial competition for host
    Rₚ::Float64          # radius of spatial competition for parasite
    ιₕ::Float64          # fitness effect of interspp interaction for host
    ιₚ::Float64          # fitness effect of interspp interaction for parasite
    Rᵢ::Float64          # radius of interspp interactions
    πₘ::Float64          # maximum probability of infection
    γ::Float64           # rate of decay of infection prb with trait difference
    αₕ::Float64          # maximum fitness effect of abiotic selection for host
    αₚ::Float64          # maximum fitness effect of abiotic selection for parasite
    Aₕ::Float64          # sensitivity of abiotic selection to trait value for host
    Aₚ::Float64          # sensitivity of abiotic selection to trait value for parasite
    θ₀ₕ::Float64         # baseline abiotic optimal trait value for host
    θ₀ₚ::Float64         # baseline abiotic optimal trait value for parasite
end

# abiotic optima as fct of location x
function θₓ(x, whch, θ₀)
    if whch == 0
        return θ₀
    end
end

# iterates through a single generation
function update(X)

    @unpack zₕ,
    zₚ,
    gₕ,
    gₚ,
    nₕ,
    nₚ,
    xₕ,
    xₚ,
    μₕ,
    μₚ,
    σₕ,
    σₚ,
    κₕ,
    κₚ,
    Rₕ,
    Rₚ,
    ιₕ,
    ιₚ,
    Rᵢ,
    πₘ,
    γ,
    αₕ,
    αₚ,
    Aₕ,
    Aₚ = X

    #
    # compute fitness for each individual
    #

    # base fitnesses for each individual
    # under absorbing boundaries
    # if individuals w/in boundaries
    # then base fitness is set to unity
    wₕ = zeros(nₕ)
    wₚ = zeros(nₚ)

    # compute distances between all individuals out here first
    htree = BallTree(transpose(xₕ))
    ptree = BallTree(transpose(xₚ))

    # accumulate effects of abiotic sel and comp on host
    for i in
        # applies only to individuals within bounds
        findall(x -> xₕ[x, 1] > 0 && xₕ[x, 2] > 0 && xₕ[x, 1] < 1 && xₕ[x, 2] < 1, 1:nₕ)

        wₕ[i] = 1

        #
        # accumulate effects of abiotic selection on host
        #
        θ = θₓ(xₕ[i], 0, θ₀ₕ)
        wₕ[i] *= αₕ * exp(-0.5 * Aₕ * (θ - zₕ[i])^2)

        #
        # accumulate effects of spatial competition on host
        #

        # # first compute vct of distances from focal ind
        # dists = zeros(nₕ)
        # for j = 1:nₕ
        #     dists[j] = sum((xₕ[i, :] .- xₕ[j, :]) .^ 2)
        # end

        # # next accumulate effects of local competition
        # wₕ[i] *= κₕ^length(findall(x -> x < Rₕ && x ≠ 0, dists))

        wₕ[i] *= κₕ^(length(inrange(htree, xₕ[i, :], Rₕ)) - 1)

        # the effect of parasitism is accounted for when we calculate parasite fitness

    end

    # accumulate effects of abiotic sel and comp on parasite
    # and effects of host-par intxns on host and parasite
    for i in
        findall(x -> xₚ[x, 1] > 0 && xₚ[x, 2] > 0 && xₚ[x, 1] < 1 && xₚ[x, 2] < 1, 1:nₚ)

        wₚ[i] = 1

        #
        # accumulate effects of abiotic selection on parasite
        #
        θ = θₓ(xₚ[i], 0, θ₀ₚ)
        wₚ[i] *= αₚ * exp(-0.5 * Aₚ * (θ - zₚ[i])^2)

        #
        # accumulate effects of spatial competition on parasite
        #

        # first compute vct of distances from focal ind
        # dists = zeros(nₚ)
        # for j = 1:nₚ
        #     dists[j] = sum((xₚ[i, :] .- xₚ[j, :]) .^ 2)
        # end

        # next accumulate effects of spatial competition
        # wₚ[i] *= κₚ^length(findall(x -> x < Rₚ && x ≠ 0, dists))
        wₚ[i] *= κₚ^(length(inrange(ptree, xₚ[i, :], Rₚ)) - 1)

        #
        # accumulate effects of host-parasite interactions
        # 

        # first compute vct of host dists from focal parasite
        # dists = zeros(nₕ)
        # for j = 1:nₕ
        #     dists[j] = sum((xₚ[i, :] .- xₕ[j, :]) .^ 2)
        # end

        # locate neighboring hosts
        nbr = inrange(htree, xₚ[i, :], Rᵢ)

        # if there's any hosts nearby
        if length(nbr) > 0

            # choose a random host within radius Rᵢ
            choice_host = rand(nbr)

            # pull trait values
            Zₕ = zₕ[choice_host]
            Zₚ = zₚ[i]

            # compute pr of infection
            πr = πₘ * exp(-γ * (Zₕ - Zₚ)^2 / 2)

            # if infection occurs, accumulate consequences
            if rand() < πr
                wₕ[choice_host] *= ιₕ
                wₚ[i] *= ιₚ
            end


        end

    end

    #
    # create next generation of hosts and then parasites
    #

    # determine number offspring for host and parasite parents
    Wₕ = zeros(nₕ)
    Wₚ = zeros(nₚ)
    for i = 1:nₕ
        Wₕ[i] = rand(Poisson(wₕ[i]))
    end
    for i = 1:nₚ
        Wₚ[i] = rand(Poisson(wₚ[i]))
    end

    # creates array of offspring
    # breeding and trait values
    # and locations

    gₕₚ = zeros(0)
    gₚₚ = zeros(0)
    x1ₕₚ = zeros(0)
    x1ₚₚ = zeros(0)
    x2ₕₚ = zeros(0)
    x2ₚₚ = zeros(0)
    for i = 1:nₕ
        append!(gₕₚ, fill(gₕ[i], Int64(Wₕ[i])))
        append!(x1ₕₚ, fill(xₕ[i, 1], Int64(Wₕ[i])))
        append!(x2ₕₚ, fill(xₕ[i, 2], Int64(Wₕ[i])))
    end
    for i = 1:nₚ
        append!(gₚₚ, fill(gₚ[i], Int64(Wₚ[i])))
        append!(x1ₚₚ, fill(xₚ[i, 1], Int64(Wₚ[i])))
        append!(x2ₚₚ, fill(xₚ[i, 2], Int64(Wₚ[i])))
    end

    # update population sizes
    nₕ = Int64(sum(Wₕ))
    nₚ = Int64(sum(Wₚ))


    # generate offspring breeding values
    gₕₚ += rand(Normal(0, √μₕ), nₕ)

    # generate offspring trait values
    zₕₚ = gₕₚ + rand(Normal(0, Eₕ), nₕ)

    # generate offspring locations (with periodic boundaries)
    # to use this, fitness calculations need to change
    # xₕₚ = mod.([x1ₕₚ x2ₕₚ] + transpose(rand( MvNormal( [0,0], Matrix(I,2,2).*σₕ ), nₕ)), 1)

    # generate offspring locations (with absorbing boundaries implemented in fitness calc above)    
    xₕₚ = [x1ₕₚ x2ₕₚ] + transpose(rand(MvNormal([0, 0], Matrix(I, 2, 2) .* σₕ), nₕ))

    gₚₚ += rand(Normal(0, √μₚ), nₚ)

    zₚₚ = gₚₚ + rand(Normal(0, Eₚ), nₚ)

    # xₚₚ = mod.([x1ₚₚ x2ₚₚ] + transpose(rand( MvNormal( [0,0], Matrix(I,2,2).*σₚ ), nₚ)), 1)    
    xₚₚ = [x1ₚₚ x2ₚₚ] + transpose(rand(MvNormal([0, 0], Matrix(I, 2, 2) .* σₚ), nₚ))

    Xₚ = hp_struct(
        zₕ = zₕₚ,
        zₚ = zₚₚ,
        gₕ = gₕₚ,
        gₚ = gₚₚ,
        nₕ = nₕ,
        nₚ = nₚ,
        xₕ = xₕₚ,
        xₚ = xₚₚ,
        μₕ = μₕ,
        μₚ = μₚ,
        σₕ = σₕ,
        σₚ = σₚ,
        Eₕ = Eₕ,
        Eₚ = Eₚ,
        θ₀ₕ = θ₀ₕ,
        θ₀ₚ = θ₀ₚ,
        κₕ = κₕ,
        κₚ = κₚ,
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

    return Xₚ

end

function sim(prs, n₀, T, hist = false)

    @unpack μₕ,
    μₚ,
    Eₕ,
    Eₚ,
    σₕ,
    σₚ,
    θ₀ₕ,
    θ₀ₚ,
    κₕ,
    κₚ,
    Rₕ,
    Rₚ,
    ιₕ,
    ιₚ,
    Rᵢ,
    πₘ,
    γ,
    αₕ,
    αₚ,
    Aₕ,
    Aₚ = prs

    nₕ = n₀[1]
    nₚ = n₀[2]

    # uniform random positions on unit square
    xₕ = rand(nₕ, 2)
    xₚ = rand(nₚ, 2)

    # initial breeding values
    gₕ = rand(Normal(θ₀ₕ, √μₕ), nₕ)
    gₚ = rand(Normal(θ₀ₚ, √μₚ), nₚ)

    # initial trait values
    if Eₕ > 0
        Eₕₘ = √Eₕ * Matrix(I, nₕ, nₕ)
        zₕ = vec(rand(MvNormal(gₕ, Eₕₘ)))
    else
        zₕ = gₕ
    end
    if Eₚ > 0
        Eₚₘ = √Eₚ * Matrix(I, nₚ, nₚ)
        zₚ = vec(rand(MvNormal(gₚ, Eₚₘ)))
    else
        zₚ = gₚ
    end

    # set up initial population
    X = hp_struct(
        zₕ = zₕ,
        zₚ = zₚ,
        gₕ = gₕ,
        gₚ = gₚ,
        nₕ = nₕ,
        nₚ = nₚ,
        xₕ = xₕ,
        xₚ = xₚ,
        μₕ = μₕ,
        μₚ = μₚ,
        Eₕ = Eₕ,
        Eₚ = Eₚ,
        σₕ = σₕ,
        σₚ = σₚ,
        θ₀ₕ = θ₀ₕ,
        θ₀ₚ = θ₀ₚ,
        κₕ = κₕ,
        κₚ = κₚ,
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

    if hist
        # set up history of population
        Xₕ = fill(X, T)

        # simulate with history
        for i = 2:T
            nh = Xₕ[i-1].nₕ
            np = Xₕ[i-1].nₚ
            if nh * np > 0
                Xₕ[i] = update(Xₕ[i-1])
            else
                Xₕ[i] = Xₕ[i-1]
            end
        end

        return Xₕ

    else

        # simulate without history
        for i = 2:T
            nh = X.nₕ
            np = X.nₚ
            if nh * np > 0
                X = update(X)            
            end
        end

        return X

    end    

end

#
# define spde pars in terms of sim pars
#

function Bₕ(γ, πₘ, ιₕ, Nₚ, Nₕ)
    # Nₚ and Nₕ need to be on same scale
    γ * πₘ * (1 - ιₕ) * Nₚ / Nₕ
end

function Bₚ(γ, πₘ, ιₚ, Nₕ)
    # Nₕ needs to be on scale of Rᵢ
    γ * πₘ * (ιₚ - 1) * (1 - exp(-Nₕ))
end

function Gₚ(μₚ, Aₚ, Bₚ)
    √(μₚ / (Aₚ + Bₚ))
end

function Gₕ(μₕ, Aₕ, Bₕ)
    √(μₕ / (Aₕ - Bₕ))
end

function v(E, G)
    E + G
end

function rₕ(αₕ, Aₕ, vₕ, πₘ, ιₕ, γ, vₚ, Nₚ, Nₕ, Bₕ)
    # Nₚ and Nₕ need to be on same scale
    log(αₕ) - 0.5 * Aₕ * vₕ +
    πₘ * (ιₕ - 1) * (1 - γ * vₚ / 2) * (Nₚ / Nₕ) +
    0.5 * Bₕ * (vₕ + vₚ)
end

function rₚ(αₚ, Aₚ, vₚ, πₘ, ιₚ, Nₕ, Bₚ, vₕ)
    # Nₕ needs to be on scale of Rᵢ
    log(αₚ) - 0.5 * Aₚ * vₚ + πₘ * (ιₚ - 1) * (1 - exp(-Nₕ)) - 0.5 * Bₚ * (vₕ + vₚ)
end

function z̄ₕ(Aₕ, θₕ, Bₕ, Aₚ, θₚ, Bₚ)
    ((Aₚ + Bₚ) * Aₕ * θₕ - Bₕ * Aₚ * θₚ) / ((Aₚ + Bₚ) * Aₕ - Bₕ * Aₚ)
end

function z̄ₚ(Aₚ, θₚ, Bₚ, Aₕ, θₕ, Bₕ)
    ((Aₕ - Bₕ) * Aₚ * θₚ + Bₚ * Aₕ * θₕ) / ((Aₕ - Bₕ) * Aₚ + Bₚ * Aₕ)
end

function mₕ(rₕ, z̄ₕ, Aₕ, θₕ, Bₕ, z̄ₚ)
    rₕ - 0.5 * (Aₕ * (θₕ - z̄ₕ)^2 + Bₕ * (z̄ₚ - z̄ₕ)^2)
end

function mₚ(rₚ, z̄ₚ, Aₚ, θₚ, Bₚ, z̄ₕ)
    rₚ - 0.5 * (Aₚ * (θₚ - z̄ₚ)^2 - Bₚ * (z̄ₕ - z̄ₚ)^2)
end

# host abundance at radius Rₕ
function Nₕ(κₕ, mₕ, μₕ, Aₕ, Bₕ)
    -(mₕ - 0.5 * √(μₕ * (Aₕ - Bₕ))) / log(κₕ)
end

# parasite abundance at radius Rₚ
function Nₚ(κₚ, mₚ, μₚ, Aₚ, Bₚ)
    -(mₚ - 0.5 * √(μₚ * (Aₚ + Bₚ))) / log(κₚ)
end

# abundance from density and radius
function Nᵨ(ρ, R)
    ρ * π * R^2
end

# density from abundance and radius
function ρ(N, R)
    N / (π * R^2)
end

# converting abundance from scale R₀ to Rₖ
function Nₖ(N₀, R₀, Rₖ)
    N₀ * (Rₖ / R₀)^2
end

function Vₚ(ρₚ, σₚ, Aₚ, Bₚ)
    1 / (σₚ^2 * ρₚ * (Aₚ + Bₚ))
end

function Vₕ(ρₕ, σₕ, Aₕ, Bₕ)
    1 / (σₕ^2 * ρₕ * (Aₕ - Bₕ))
end

function ξₚ(Gₚ, σₚ, Aₚ, Bₚ)
    σₚ / √(Gₚ * (Aₚ + Bₚ))
end

function ξₕ(Gₕ, σₕ, Aₕ, Bₕ)
    σₕ / √(Gₕ * (Aₕ - Bₕ))
end

