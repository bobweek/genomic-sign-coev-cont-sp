################################################################################
##
## AUTHOR: Bob Week
##
## DATE: 08/23/2021
##
## In this script we provide definitions of data structures for model parameters
## and for state variables used in our simulations of individual-based models.
##
################################################################################

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

# abiotic optima as fct of location x
function θₓ(x,whch,θ₀)
    if whch==0
        return θ₀
    end
end

# iterates through a single generation
function update(X)

    @unpack zₕ, zₚ, gₕ, gₚ, nₕ, nₚ, xₕ, xₚ, μₕ, μₚ, σₕ, σₚ, 
        κₕ, κₚ, Rₕ, Rₚ, ιₕ, ιₚ, Rᵢ, πₘ, γ, αₕ, αₚ, Aₕ, Aₚ = X

    #
    # compute fitness for each individual
    #

    # base fitnesses for each individual
    wₕ = ones(nₕ)
    wₚ = ones(nₚ)

    # accumulate effects of abiotic sel and comp on host
    for i in 1:nₕ

        #
        # accumulate effects of abiotic selection on host
        #
        θ = θₓ(xₕ[i],0,θ₀ₕ)
        wₕ[i] *= αₕ*exp(-Aₕ*(θ-zₕ[i])^2)

        #
        # accumulate effects of spatial competition on host
        #

        # first compute vct of distances from focal ind
        dists = zeros(nₕ)
        for j in 1:nₕ
            dists[j] = sum((xₕ[i].-xₕ[j]).^2)
        end

        # next accumulate effects of spatial competition
        wₕ[i] *= κₕ^length(findall(x -> x<Rₕ && x≠0, dists))
        
    end

    # accumulate effects of abiotic sel and comp on parasite
    # and effects of host-par intxns on host and parasite
    for i in 1:nₚ

        #
        # accumulate effects of abiotic selection on parasite
        #
        θ = θₓ(xₚ[i],0,θ₀ₚ)
        wₚ[i] *= αₚ*exp(-Aₚ*(θ-zₚ[i])^2)

        #
        # accumulate effects of spatial competition on parasite
        #

        # first compute vct of distances from focal ind
        dists = zeros(nₚ)
        for j in 1:nₚ
            dists[j] = sum((xₚ[i].-xₚ[j]).^2)
        end

        # next accumulate effects of spatial competition
        wₚ[i] *= κₚ^length(findall(x -> x<Rₚ && x≠0, dists))

        #
        # accumulate effects of host-parasite interactions
        # 

        # first compute vct of host dists from focal parasite
        dists = zeros(nₕ)
        for j in 1:nₕ
            dists[j] = sum((xₚ[i].-xₕ[j]).^2)
        end

        # next choose a random host within radius Rᵢ
        choice_host = rand(findall(x -> x<Rᵢ, dists))

        # pull trait values
        Zₕ = zₕ[choice_host]
        Zₚ = zₚ[i]

        # compute pr of infection
        π = πₘ*exp(-γ*(Zₕ-Zₚ)^2)

        # if infection occurs, accumulate consequences
        if(rand()<π)
            wₕ[choice_host] *= ιₕ
            wₚ[i] *= ιₚ
        end
        
    end

    #
    # create next generation of hosts and then parasites
    #

    # determine number offspring for host and parasite parents
    Wₕ = zeros(nₕ)
    Wₚ = zeros(nₚ)
    for i in 1:nₕ
        Wₕ[i] = rand(Poisson(wₕ[i]))
    end
    for i in 1:nₚ
        Wₚ[i] = rand(Poisson(wₚ[i]))
    end

    # creates array of offspring
    # breeding and trait values
    # and location

    gₕₚ = zeros(0)
    gₚₚ = zeros(0)    
    x1ₕₚ = zeros(0)
    x1ₚₚ = zeros(0)
    x2ₕₚ = zeros(0)
    x2ₚₚ = zeros(0)    
    for i in 1:nₕ
        append!(gₕₚ,fill(gₕ[i],Int64(Wₕ[i])))
        append!(x1ₕₚ,fill(xₕ[i,1],Int64(Wₕ[i])))
        append!(x2ₕₚ,fill(xₕ[i,2],Int64(Wₕ[i])))
    end
    for i in 1:nₚ
        append!(gₚₚ,fill(gₚ[i],Int64(Wₚ[i])))
        append!(x1ₚₚ,fill(xₚ[i,1],Int64(Wₚ[i])))
        append!(x2ₚₚ,fill(xₚ[i,2],Int64(Wₚ[i])))
    end

    # update population sizes
    nₕ = Int64(sum(Wₕ))
    nₚ = Int64(sum(Wₚ))

    if nₕ>0
        
        # generate offspring breeding values
        gₕₚ += rand( Normal( 0, √μₕ ), nₕ)

        # generate offspring trait values
        zₕₚ = gₕₚ + rand( Normal( 0, Eₕ ), nₕ)
    
        # generate offspring locations (with periodic boundaries)
        xₕₚ = mod.([x1ₕₚ x2ₕₚ] + transpose(rand( MvNormal( [0,0], Matrix(I,2,2).*σₕ ), nₕ)), 1)
    
    end

    if nₚ>0
    
        gₚₚ += rand( Normal( 0, √μₚ ), nₚ)

        zₚₚ = gₚₚ + rand( Normal( 0, Eₚ ), nₚ)

        xₚₚ = mod.([x1ₚₚ x2ₚₚ] + transpose(rand( MvNormal( [0,0], Matrix(I,2,2).*σₚ ), nₚ)), 1)

    end
    
    Xₚ = hp_struct(zₕ=zₕₚ, zₚ=zₚₚ, gₕ=gₕₚ, gₚ=gₚₚ, nₕ=nₕ, nₚ=nₚ, xₕ=xₕₚ, xₚ=xₚₚ, μₕ=μₕ, μₚ=μₚ, σₕ=σₕ, σₚ=σₚ, Eₕ=Eₕ, Eₚ=Eₚ,
        θ₀ₕ=θ₀ₕ, θ₀ₚ=θ₀ₚ, κₕ=κₕ, κₚ=κₚ, Rₕ=Rₕ, Rₚ=Rₚ, ιₕ=ιₕ, ιₚ=ιₚ, Rᵢ=Rᵢ, πₘ=πₘ, γ=γ, αₕ=αₕ, αₚ=αₚ, Aₕ=Aₕ, Aₚ=Aₚ)

    return Xₚ

end
