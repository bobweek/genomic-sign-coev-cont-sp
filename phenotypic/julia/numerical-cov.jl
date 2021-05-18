using SpecialFunctions     # provides Bessel and Gamma functions
using GaussianRandomFields # for simulating (univariate) random fields
using DSP                  # for numerical convolution
using Plots
using Parameters
pyplot()

#
# to visualize host-parasite cross-covariance, we need to numerically convolve Bessel functions
#
# in particular, we want Cₕₚ(x) = [K₀(⋅)∗M(⋅|ν,κ)](x)
#

# data type that holds model parameters
@with_kw mutable struct CoevPars

    # primary parameters
    Gₕ::Float64 # additive genetic variance of host
	Gₚ::Float64 # additive genetic variance of parasite
	Nₕ::Float64 # effective population size of host
    Nₚ::Float64 # effective population size of parasite
	Aₕ::Float64	# strength of abiotic selection on host
	Aₚ::Float64	# strength of abiotic selection on parasite
	Bₕ::Float64	# strength of biotic selection on host
	Bₚ::Float64	# strength of biotic selection on parasite
	Dₕ::Float64	# host disperal distance
	Dₚ::Float64	# parasite disperal distance

    # parameters for calculating local adaptation measures
    σₕ²::Float64 # expressed variance of host
	σₚ²::Float64 # expressed variance of parasite
	rₕ::Float64  # intrinsic growth rate of host
    rₚ::Float64  # intrinsic growth rate of parasite
	
end

# covariance functions
Cₕₕ = function (x,p)

    # x = distance
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,Dₕ,Dₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = √(Dₕ / (2 * Gₕ * (Aₕ - Bₕ)))
    ξₚ = √(Dₚ / (2 * Gₚ * (Aₚ + Bₚ)))

    # variance in mean trait values
    # due to random genetic drift
    local_varₕ = 1 / (Nₕ * Dₕ * (Aₕ - Bₕ))
    local_varₚ = 1 / (Nₚ * Dₚ * (Aₚ + Bₚ))

    return(local_varₕ * x * besselk(1, x / ξₕ) / ξₕ)
end

Cₚₚ = function (x,p)

    # x = distance
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,Dₕ,Dₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = √(Dₕ / (2 * Gₕ * (Aₕ - Bₕ)))
    ξₚ = √(Dₚ / (2 * Gₚ * (Aₚ + Bₚ)))

    # variance in mean trait values
    # due to random genetic drift
    local_varₕ = 1 / (Nₕ * Dₕ * (Aₕ - Bₕ))
    local_varₚ = 1 / (Nₚ * Dₚ * (Aₚ + Bₚ))

    return(local_varₚ * x * besselk(1, x / ξₚ) / ξₚ)
end

# now for the convolution...
Cₕₚ = function (l₁, u₁, s₁, l₂, u₂, s₂, p)

    # lᵢ is lower bound for i-th coordinate
    # uᵢ is upper bound for i-th coordinate
    # sᵢ is step size for i-th coordinate
    # p = model parameters

    # define range to convolve over
    x₁ = l₁:s₁:u₁
    x₂ = l₂:s₂:u₂
    x₁ = filter(x -> x ≠ 0, x₁) # avoids singularities
    x₂ = filter(x -> x ≠ 0, x₂)

    # unpack model parameters
    @unpack Gₕ,Gₚ,Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,Dₕ,Dₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = √(Dₕ / (2 * Gₕ * (Aₕ - Bₕ)))
    ξₚ = √(Dₚ / (2 * Gₚ * (Aₚ + Bₚ)))

    # variance in mean trait values
    # due to random genetic drift
    local_varₕ = 1 / (Nₕ * Dₕ * (Aₕ - Bₕ))
    local_varₚ = 1 / (Nₚ * Dₚ * (Aₚ + Bₚ))


    # holds values of Bessel K₀ function with host parameters
    K₀ₕ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            K₀ₕ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = besselk(0, u / ξₕ)
        end
    end

    # holds values of Bessel K₀ function with parasite parameters
    K₀ₚ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            K₀ₚ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = besselk(0, u / ξₚ)
        end
    end

    M₁ₕ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            M₁ₕ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = u * besselk(1, u / ξₕ) / ξₕ
        end
    end

    M₁ₚ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            M₁ₚ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = u * besselk(1, u / ξₚ) / ξₚ
        end
    end

    K₀ₚM₁ₕ = conv(K₀ₚ, M₁ₕ)
    K₀ₕM₁ₚ = conv(K₀ₕ, M₁ₚ)

    Cₕₚ = (2 .* Gₚ .* Bₚ .* local_varₕ .* (K₀ₚM₁ₕ) ./ Dₚ) .- (2 .* Gₕ .* Bₕ .* local_varₚ .* (K₀ₕM₁ₚ) ./ Dₕ)

    return(Cₕₚ)

end

# plots spatial correlations and cross-correlation
plotSpCorr = function (m,s,p)

    # m = max distance
    # s = step size (ie., resolution)
    # p = model parameters

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m)

    CHH = zeros(length(U))
    CPP = zeros(length(U))
    count = 1
    for u in U
        CHH[count] = Cₕₕ(u,p)
        CPP[count] = Cₚₚ(u,p)
        count += 1
    end

    # unpack some parameters
    @unpack Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,Dₕ,Dₚ = p

    # calculate local variances to normalize correlations
    local_varₕ = 1 / (Nₕ * Dₕ * (Aₕ - Bₕ))
    local_varₚ = 1 / (Nₚ * Dₚ * (Aₚ + Bₚ))

    ttle = string("Dispersal Ratio: σₕ/σₚ = ",√(Dₕ/Dₚ))

    # make array with correlations on columns
    L = length(U)-1
    ρs = hcat( CHH, CPP, CHP[(L-1):(2 * L - 1),L]/sqrt(local_varₕ*local_varₚ) )
    plot(U, ρs,label=["ρₕ" "ρₚ" "ρₕₚ"],title=ttle, legendfontsize=10)
    ylabel!("Trait Correlation")
    xlabel!("Spatial Lag")

end

# plots measures of local adaptation
plotLocAdapt = function (m,s,p,type)

    # m = max distance
    # s = step size (ie., resolution)
    # p = model parameters

    # unpack some parameters
    @unpack Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,Dₕ,Dₚ,rₕ,rₚ,σₕ²,σₚ² = p
    
    # these are same as Cₕₕ(0,p) and Cₚₚ(0,p) resp.
    # local_varₕ = 1 / (Nₕ * Dₕ * (Aₕ - Bₕ))
    # local_varₚ = 1 / (Nₚ * Dₚ * (Aₚ + Bₚ))

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m)

    L = length(U)-1
    CHP0 = CHP[L,L]
    CHPX = CHP[(L-1):(2 * L - 1),L]

    Δₕ = Bₕ .* (CHPX .- CHP0)
    Δₚ = Bₚ .* (CHP0 .- CHPX)
    # Δₕₚ = (rₕ-rₚ) .+ (Bₚ-Bₕ).*CHPX .+ ( (Bₕ-Aₕ-Bₚ)*σₕ² + (Aₚ+Bₚ+Bₕ)*σₚ² + (Bₚ-Aₕ-Bₕ)*local_varₕ + (Aₚ+Bₚ-Bₕ)*local_varₚ )/2    

    if type == "D"
        ttle = string("Dispersal Ratio: σₕ/σₚ = ",√(Dₕ/Dₚ))
    elseif type == "B"
        ttle = string("Coevolution Ratio: Bₕ/Bₚ = ",Bₕ/Bₚ)
    end

    # make array with LA measures on columns    
    las = hcat(Δₕ, Δₚ)
    plot(U, las,label=["H" "P"],title=ttle, legendfontsize=10)
    ylabel!("Fitness Difference")
    xlabel!("Spatial Lag")

end

# same dispersal distances
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, Dₕ = 10^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
eqLA = plotLocAdapt(60,0.1,p,"D")
title!("")
# plot(eqCorr, eqLA, xlims=(0,40), layout = (2,1), size = (400,600))
plot(eqLA, xlims=(0,40), title = "", size = (500,400))
savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/eq-pars.png")

# parasite disperses further than host
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, Dₕ = 10^2, Dₚ = 100^2, rₕ = 0, rₚ = 0)
# Corrp = plotSpCorr(60,0.1,p)
LAp = plotLocAdapt(60,0.1,p,"D")
# plot(Corrp, LAp, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/p-further.png")


# host disperses a tiny bit further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, Dₕ = 10.5^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# Corrhl = plotSpCorr(60,0.1,p)
LAht = plotLocAdapt(60,0.1,p,"D")
# plot(Corrhl, LAhl, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-l-further.png")

# host disperses a little bit further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, Dₕ = 20^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# Corrhl = plotSpCorr(60,0.1,p)
LAhl = plotLocAdapt(60,0.1,p,"D")
# plot(Corrhl, LAhl, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-l-further.png")

# host disperses an intermediate bit further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, Dₕ = 50^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# Corrhi = plotSpCorr(60,0.1,p)
LAhi = plotLocAdapt(60,0.1,p,"D")
# plot(Corrhi, LAhi, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-i-further.png")

# host disperses much further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, Dₕ = 100^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# Corrh = plotSpCorr(60,0.1,p)
LAh = plotLocAdapt(60,0.1,p,"D")
# plot(Corrh, LAh, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-further.png")

# plot(Corrp, eqCorr, Corrhl, Corrhi, Corrh, eqLA, LAp, LAhl, LAhi, LAh, xlims=(0,40), layout=(2,5), size=(1600,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/corrs-las.png")

# equal strength
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, Dₕ = 10^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
eqBLA = plotLocAdapt(60,0.1,p,"B")

# host stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.02, Bₚ = 0.01, Dₕ = 10^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
hsLA = plotLocAdapt(60,0.1,p,"B")

# host little stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.0106, Bₚ = 0.01, Dₕ = 10^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
hlsLA = plotLocAdapt(60,0.1,p,"B")

# host much stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.1, Bₚ = 0.01, Dₕ = 10^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
hmsLA = plotLocAdapt(60,0.1,p,"B")

# parasite stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.02, Dₕ = 10^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
psLA = plotLocAdapt(60,0.1,p,"B")

# parasite much stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, σₕ²=10, σₚ²=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.1, Dₕ = 10^2, Dₚ = 10^2, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
pmsLA = plotLocAdapt(60,0.1,p,"B")

plot(LAp, eqLA, LAht, LAhl, LAhi, LAh, xlims=(0,40), size=(1200,600))
savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/d-las.png")

plot(pmsLA, eqBLA, hlsLA, hmsLA, xlims=(0,40), size=(1200,600))
savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/b-las.png")

#
# should change strengths of biotic selection to see how
# that changes spatial correlations
#
