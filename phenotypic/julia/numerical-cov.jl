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
	σₕ::Float64	# host disperal distance
	σₚ::Float64	# parasite disperal distance

    # parameters for calculating local adaptation measures
    vₕ::Float64 # expressed variance of host
	vₚ::Float64 # expressed variance of parasite
	rₕ::Float64  # intrinsic growth rate of host
    rₚ::Float64  # intrinsic growth rate of parasite
	
end

# covariance functions
Cₕₕ = function (x,p)

    # x = distance
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = σₕ / √(Gₕ * (Aₕ - Bₕ))

    # variance in mean trait values
    Vₕ = 1 / (Nₕ * σₕ * (Aₕ - Bₕ))

    return(local_varₕ * x * besselk(1,√2 * x / ξₕ) / ξₕ)
end

Cₚₚ = function (x,p)

    # x = distance
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₚ = σₚ / √(Gₚ * (Aₚ + Bₚ))

    # variance in mean trait values
    Vₚ = 1 / (Nₚ * σₚ * (Aₚ + Bₚ))

    return(Vₚ * x * besselk(1,√2 * x / ξₚ) / ξₚ)
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
    @unpack Gₕ,Gₚ,Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = σₕ / √(Gₕ * (Aₕ - Bₕ))
    ξₚ = σₚ / √(Gₚ * (Aₚ + Bₚ))

    # variance in mean trait values
    Vₚ = 1 / (Nₚ * σₚ * (Aₚ + Bₚ))
    Vₕ = 1 / (Nₕ * σₕ * (Aₕ - Bₕ))

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
            M₁ₕ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = √2 * u * besselk(1,√2 * u / ξₕ) / ξₕ
        end
    end

    M₁ₚ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            M₁ₚ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = √2 * u * besselk(1,√2 * u / ξₚ) / ξₚ
        end
    end

    K₀ₚM₁ₕ = conv(K₀ₚ, M₁ₕ)
    K₀ₕM₁ₚ = conv(K₀ₕ, M₁ₚ)

    Cₕₚ = (2 .* Gₚ .* Bₚ .* Vₕ .* (K₀ₚM₁ₕ) ./ σₚ) .- (2 .* Gₕ .* Bₕ .* Vₚ .* (K₀ₕM₁ₚ) ./ σₕ)

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
    @unpack Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # calculate local variances to normalize correlations
    local_varₕ = 1 / (Nₕ * σₕ * (Aₕ - Bₕ))
    local_varₚ = 1 / (Nₚ * σₚ * (Aₚ + Bₚ))

    ttle = string("Dispersal Ratio: σₕ/σₚ = ", σₕ/σₚ)

    # make array with correlations on columns
    L = length(U)-1
    ρs = hcat( CHH, CPP, CHP[(L-1):(2 * L - 1),L]/sqrt(Vₕ*Vₚ) )
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
    @unpack Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ,rₕ,rₚ,vₕ,vₚ = p

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m)

    L = length(U)-1
    CHP0 = CHP[L,L]
    CHPX = CHP[(L-1):(2 * L - 1),L]

    Δₕ = Bₕ .* (CHPX .- CHP0)
    Δₚ = Bₚ .* (CHP0 .- CHPX)

    if type == "σ"
        ttle = string("Dispersal Ratio: σₕ/σₚ = ",σₕ/σₚ)
    elseif type == "B"
        ttle = string("Coevolution Ratio: Bₕ/Bₚ = ",Bₕ/Bₚ)
    end

    # make array with LA measures on columns    
    las = hcat(Δₕ, Δₚ)
    plot(U, las,label=["H" "P"],title=ttle, legendfontsize=10)
    ylabel!("Fitness Difference")
    xlabel!("Spatial Lag")

end


# returns a measure of local adaptation considering limited dispersal
LimDispLA = function (m,s,p,type)

    # m = max distance
    # s = step size (ie., resolution)
    # p = model parameters

    # unpack some parameters
    @unpack Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ,rₕ,rₚ,vₕ,vₚ = p
    
    # these are same as Cₕₕ(0,p) and Cₚₚ(0,p) resp.
    # local_varₕ = 1 / (Nₕ * σₕ * (Aₕ - Bₕ))
    # local_varₚ = 1 / (Nₚ * σₚ * (Aₚ + Bₚ))

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m)

    L = length(U)-1
    CHP0 = CHP[L,L]
    CHPX = CHP[(L-1):(2 * L - 1),L] # replace this one with the value at the expected distance

    Δₕ = Bₕ .* (CHPX .- CHP0)
    Δₚ = Bₚ .* (CHP0 .- CHPX)
    # Δₕₚ = (rₕ-rₚ) .+ (Bₚ-Bₕ).*CHPX .+ ( (Bₕ-Aₕ-Bₚ)*vₕ + (Aₚ+Bₚ+Bₕ)*vₚ + (Bₚ-Aₕ-Bₕ)*Vₕ + (Aₚ+Bₚ-Bₕ)*Vₚ )/2    

    if type == "D"
        ttle = string("Dispersal Ratio: σₕ/σₚ = ",σₕ/σₚ)
    elseif type == "B"
        ttle = string("Coevolution Ratio: Bₕ/Bₚ = ",Bₕ/Bₚ)
    end

    # make array with LA measures on columns    
    las = hcat(Δₕ, Δₚ)
    plot(U, las,label=["H" "P"],title=ttle, legendfontsize=10)
    ylabel!("Fitness Difference")
    xlabel!("Spatial Lag")

end

# returns a measure of local adaptation considering limited dispersal
ClassicLA = function (m,s,p)

    # m = max distance
    # s = step size (ie., resolution)
    # p = model parameters

    # unpack some parameters
    @unpack Nₕ,Nₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ,rₕ,rₚ,vₕ,vₚ = p
    
    # these are same as Cₕₕ(0,p) and Cₚₚ(0,p) resp.
    # local_varₕ = 1 / (Nₕ * σₕ * (Aₕ - Bₕ))
    # local_varₚ = 1 / (Nₚ * σₚ * (Aₚ + Bₚ))

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m)

    L = length(U)-1
    CHP0 = CHP[L,L]
    CHPX = CHP[(L-1):(2 * L - 1),L] # set to zero

    
    ℒₕ = -Bₕ .* CHP0
    ℒₚ =  Bₚ .* CHP0

    ℒ = [ℒₕ, ℒₚ]

    return(ℒ)

end