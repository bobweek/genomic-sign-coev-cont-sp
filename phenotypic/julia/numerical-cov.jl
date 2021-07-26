using SpecialFunctions     # provides Bessel and Gamma functions
using DSP                  # for numerical convolution
using Cuba
using Plots
using Parameters
gr()

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
	ρₕ::Float64 # effective population size of host
    ρₚ::Float64 # effective population size of parasite
	Aₕ::Float64	# strength of abiotic selection on host
	Aₚ::Float64	# strength of abiotic selection on parasite
	Bₕ::Float64	# strength of biotic selection on host
	Bₚ::Float64	# strength of biotic selection on parasite
	σₕ::Float64	# host disperal distance
	σₚ::Float64	# parasite disperal distance

    # other parameters
    vₕ::Float64  # expressed variance of host
	vₚ::Float64  # expressed variance of parasite
	rₕ::Float64  # intrinsic growth rate of host
    rₚ::Float64  # intrinsic growth rate of parasite
    𝓝ₕ::Float64 # Wright's neighborhood size of host
    𝓝ₚ::Float64 # Wright's neighborhood size of parasite
	
end

# power spectra
Sₕ = function(k,p)

    if length(k) != 2
        print("k needs to be a 2d vector")
        return
    end

    # k = [k₁, k₂] = wavenumber
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    knorm = sqrt(k[1]^2+k[2]^2)
    num = Bₕ^2*Gₕ^2*Gₚ/ρₚ + Gₕ*( Gₚ*(Aₚ+Bₚ)+0.5*σₚ^2*knorm^2 )^2/ρₕ
    den = (2*π)*( Bₕ*Bₚ*Gₕ*Gₚ + (Gₕ*(Aₕ-Bₕ)+0.5*(σₕ*knorm)^2) * (Gₚ*(Aₚ+Bₚ)+0.5*(σₚ*knorm)^2) )^2

    return(num/den)

end

Sₚ = function(k,p)

    if length(k) != 2
        print("k needs to be a 2d vector")
        return
    end

    # k = [k₁, k₂] = wavenumber
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    knorm = sqrt(k[1]^2+k[2]^2)
    num = Bₚ^2*Gₚ^2*Gₕ/ρₕ + Gₚ*( Gₕ*(Aₕ-Bₕ)+0.5*σₕ^2*knorm^2 )^2/ρₚ
    den = (2*π)*( Bₕ*Bₚ*Gₕ*Gₚ + (Gₕ*(Aₕ-Bₕ)+0.5*(σₕ*knorm)^2) * (Gₚ*(Aₚ+Bₚ)+0.5*(σₚ*knorm)^2) )^2

    return(num/den)

end

Sₕₚ = function(k,p)

    if length(k) != 2
        print("k needs to be a 2d vector")
        return
    end

    # k = [k₁, k₂] = wavenumber
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    knorm = sqrt(k[1]^2+k[2]^2)
    num = Bₚ*Gₚ*Gₕ*( Gₚ*(Aₚ+Bₚ)+0.5*(σₚ*knorm)^2 )/ρₕ - Bₕ*Gₕ*Gₚ*( Gₕ*(Aₕ-Bₕ)+0.5*(σₕ*knorm)^2 )/ρₚ
    den = ( Bₕ*Bₚ*Gₕ*Gₚ + (Gₕ*(Aₕ-Bₕ)+0.5*(σₕ*knorm)^2) * (Gₚ*(Aₚ+Bₚ)+0.5*(σₚ*knorm)^2) )^2

    return(num/den)

end

S̃ₕₚ = function(k,p)

    if length(k) != 2
        print("k needs to be a 2d vector")
        return
    end

    # k = [k₁, k₂] = wavenumber
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p
    
    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = σₕ / √(Gₕ * (Aₕ - Bₕ))
    ξₚ = σₚ / √(Gₚ * (Aₚ + Bₚ))

    # variance in mean trait values
    Vₕ = 1 / (ρₕ * σₕ^2 * (Aₕ - Bₕ))
    Vₚ = 1 / (ρₚ * σₚ^2 * (Aₚ + Bₚ))
    
    knorm = sqrt(k[1]^2+k[2]^2)

    numₚ = Bₚ*Vₕ*ξₕ^2
    deρₚ = (Aₚ+Bₚ)*(1+0.5*(ξₕ*knorm)^2)^2*(1+0.5*(ξₚ*knorm)^2)

    numₕ = Bₕ*Vₚ*ξₚ^2
    deρₕ = (Aₕ-Bₕ)*(1+0.5*(ξₚ*knorm)^2)^2*(1+0.5*(ξₕ*knorm)^2)

    return((numₚ/deρₚ) - (numₕ/deρₕ))

end

CubaSₕₚ = function(y,p)

    if length(y) != 2
        print("y needs to be a 2d vector")
        return
    end

    k₁  = (2*y[1]-1) / (y[1]*(1-y[1]))
    dk₁ = (2*y[1]^2-2*y[1]+1) / (y[1]^2*(1-y[1])^2)

    k₂  = (2*y[2]-1) / (y[2]*(1-y[2]))
    dk₂ = (2*y[2]^2-2*y[2]+1) / (y[2]^2*(1-y[2])^2)
    
    # k = [k₁, k₂] = wavenumber
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    knorm = sqrt(k₁^2+k₂^2)
    num₁ = Bₚ*Gₚ*Gₕ*( Gₚ*(Aₚ+Bₚ)+0.5*(σₚ*knorm)^2 )/ρₕ
    num₂ = Bₕ*Gₕ*Gₚ*( Gₕ*(Aₕ-Bₕ)+0.5*(σₕ*knorm)^2 )/ρₚ
    num = num₁ - num₂
    den = ( Bₕ*Bₚ*Gₕ*Gₚ + (Gₕ*(Aₕ-Bₕ)+0.5*(σₕ*knorm)^2) * (Gₚ*(Aₚ+Bₚ)+0.5*(σₚ*knorm)^2) )^2

    return(dk₁*dk₂*num/den)

end

CubaD̂ = function(y,p)
    
    if length(y) != 2
        print("y needs to be a 2d vector")
        return
    end

    k₁  = (2*y[1]-1) / (y[1]*(1-y[1]))
    dk₁ = (2*y[1]^2-2*y[1]+1) / (y[1]^2*(1-y[1])^2)

    k₂  = (2*y[2]-1) / (y[2]*(1-y[2]))
    dk₂ = (2*y[2]^2-2*y[2]+1) / (y[2]^2*(1-y[2])^2)
    
    D̂ = exp(-0.5*(k₁^2+k₂^2)*(p.σₕ^2+p.σₚ^2))/(2*π)

    return(D̂*dk₁*dk₂)

end

# calculating expected covariance via numerical integration wrt distance distribution in fourier space
C̄ₕₚ = function (p)

    cout = cuhre((k,f) -> f[1] = CubaD̂(k,p)*CubaSₕₚ(k,p)/(2*π), 2)

    if cout[5]!=0
        print("error\n")
        return(0)
    end

    C̄ = cout[1][1]

    return(C̄)
    
end

# covariance functions
Cₕₕ = function (x,p)

    # x = distance
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = σₕ / √(Gₕ * (Aₕ - Bₕ))

    # variance in mean trait values
    Vₕ = 1 / (ρₕ * σₕ^2 * (Aₕ - Bₕ))

    return(Vₕ * √2 * x * besselk(1,√2 * x / ξₕ) / ξₕ)
end

Cₚₚ = function (x,p)

    # x = distance
    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₚ = σₚ / √(Gₚ * (Aₚ + Bₚ))

    # variance in mean trait values
    Vₚ = 1 / (ρₚ * σₚ^2 * (Aₚ + Bₚ))

    return(Vₚ * √2 * x * besselk(1,√2 * x / ξₚ) / ξₚ)
end

# the marginal covariance
Cₕₚ₀ = function (p)

    # p = model parameters

    # unpack model parameters
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = σₕ / √(Gₕ * (Aₕ - Bₕ))
    ξₚ = σₚ / √(Gₚ * (Aₚ + Bₚ))

    numₚ = Bₚ*Gₕ*Gₚ*(ξₕ^6*ξₚ^2-ξₕ^4*ξₚ^4-log(ξₕ^2)*ξₕ^4*ξₚ^4+log(ξₚ^2)*ξₕ^4*ξₚ^4)
    deρₚ = ρₕ*(ξₕ^2-ξₚ^2)^2*σₕ^4*σₚ^2

    numₕ = Bₕ*Gₕ*Gₚ*(ξₕ^4*ξₚ^4-log(ξₕ^2)*ξₕ^4*ξₚ^4+log(ξₚ^2)*ξₕ^4*ξₚ^4-ξₕ^2*ξₚ^6)
    deρₕ = ρₕ*(ξₕ^2-ξₚ^2)^2*σₕ^4*σₚ^2

    # CHP0 = 8*Gₕ*Gₚ*(ξₕ*ξₚ)^2 * ( Bₚ*(ξₕ^4+(ξₕ*ξₚ)^2*(2*log(ξₚ/ξₕ)-1))/(ρₕ*σₕ^2) - Bₕ*(ξₚ^4+(ξₕ*ξₚ)^2*(2*log(ξₕ/ξₚ)-1))/(ρₚ*σₚ^2) ) / ( σₕ^2*σₚ^2*(ξₕ^2-ξₚ^2)^2 )

    return((numₚ/deρₚ)+(numₕ/deρₕ))

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
    @unpack Gₕ,Gₚ,ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # compute characteristic lengths of intraspecific spatial variation
    ξₕ = σₕ / √(Gₕ * (Aₕ - Bₕ))
    ξₚ = σₚ / √(Gₚ * (Aₚ + Bₚ))

    # variance in mean trait values
    Vₚ = 1 / (ρₚ * σₚ^2 * (Aₚ + Bₚ))
    Vₕ = 1 / (ρₕ * σₕ^2 * (Aₕ - Bₕ))

    # holds values of Bessel K₀ function with host parameters
    K₀ₕ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            K₀ₕ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = 2 .* besselk(0, √2*u/ξₕ) / ξₕ^2
        end
    end

    # holds values of Bessel K₀ function with parasite parameters
    K₀ₚ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            K₀ₚ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = 2 .* besselk(0, √2*u/ξₚ) / ξₚ^2
        end
    end

    M₁ₕ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            M₁ₕ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = Vₕ * √2 * u * besselk(1, √2*u/ξₕ) / ξₕ
        end
    end

    M₁ₚ = zeros(length(x₁), length(x₂))
    for x in x₁
        for y in x₂
            u = √(x^2 + y^2)
            M₁ₚ[findfirst(z -> z == x, x₁),findfirst(z -> z == y, x₂)] = Vₚ * √2 * u * besselk(1, √2*u/ξₚ) / ξₚ
        end
    end

    K₀ₚM₁ₕ = conv(K₀ₚ, M₁ₕ)
    K₀ₕM₁ₚ = conv(K₀ₕ, M₁ₚ)


    CHP0 = Cₕₚ₀(p)

    Cₕₚ = ( Bₚ.*(K₀ₚM₁ₕ)./(Aₚ+Bₚ) ) .- ( Bₕ.*(K₀ₕM₁ₚ)./(Aₕ-Bₕ) )

    Cₕₚ *= CHP0/maximum(Cₕₚ)

    return(Cₕₚ)

end

# the marginal covariance computed numerically
Cₕₚ₀NUM = function (m,s,p)

    # m = max distance
    # s = step size (ie., resolution)
    # p = model parameters

    # unpack some parameters
    @unpack ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ,rₕ,rₚ,vₕ,vₚ = p
    
    # these are same as Cₕₕ(0,p) and Cₚₚ(0,p) resp.
    # local_varₕ = 1 / (ρₕ * σₕ * (Aₕ - Bₕ))
    # local_varₚ = 1 / (ρₚ * σₚ * (Aₚ + Bₚ))

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m)

    L = length(U)-1
    CHP0 = CHP[L,L]
    
    return(CHP0)

end

# the covariance at expected distance
Cₕₚd̄ = function (m,s,p)

    # m = max distance
    # s = step size (ie., resolution)
    # p = model parameters

    # unpack some parameters
    @unpack ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ,rₕ,rₚ,vₕ,vₚ = p
    
    # expected distance
    d̄ = √(π*(σₕ^2+σₚ^2)/2)

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m)

    ind_d̄ = findfirst( z -> z >= d̄, U)

    L = length(U) - 1
    Ld̄ = ind_d̄ + length(U) - 1
    CHPd̄ = CHP[Ld̄,L]
    
    return(CHPd̄)

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
    @unpack ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ = p

    # calculate local variances to normalize correlations
    Vₕ = 1 / (ρₕ * σₕ^2 * (Aₕ - Bₕ))
    Vₚ = 1 / (ρₚ * σₚ^2 * (Aₚ + Bₚ))

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
    @unpack ρₕ,ρₚ,Aₕ,Aₚ,Bₕ,Bₚ,σₕ,σₚ,rₕ,rₚ,vₕ,vₚ = p

    # compute cross-covariance function
    CHP = Cₕₚ(-m,m,s,-m,m,s,p)

    # define range
    U = 0:s:(2*m-s)

    # expected distance
    # integrate this into the plot
    d̄ = √(π*(σₕ^2+σₚ^2)/2)

    L = length(U)
    CHP0 = CHP[L,L]
    CHPX = CHP[L:(2*L-1),L]

    ℓₕ = Bₕ .* (CHPX .- CHP0)
    ℓₚ = Bₚ .* (CHP0 .- CHPX)

    if type == "σ"
        ttle = string("Dispersal Ratio: σₕ/σₚ = ",σₕ/σₚ)
    elseif type == "B"
        ttle = string("Coevolution Ratio: Bₕ/Bₚ = ",Bₕ/Bₚ)
    end

    # make array with LA measures on columns    
    las = hcat(ℓₕ, ℓₚ)
    plot(U, las,label=["H" "P"],title=ttle, legendfontsize=10)
    ylabel!("Fitness Difference")
    xlabel!("Spatial Lag")

end


# returns a measure of local adaptation that accounts for limited dispersal
# slow and only approximate, replacing with below
# ℓₗᵢₘ = function (m,s,p)

#     # m = max distance
#     # s = step size (ie., resolution)
#     # p = model parameters

#     # unpack some parameters
#     @unpack Bₕ,Bₚ = p

#     # compute cross-covariance function
#     CHPd̄ = Cₕₚd̄(m,s,p)
#     CHP0 = Cₕₚ₀(p)

#     ℓₕ = Bₕ * (CHPd̄ - CHP0)
#     ℓₚ = Bₚ * (CHP0 - CHPd̄)
    
#     ℓ = [ℓₕ, ℓₚ]

#     return(ℓ)

# end


# returns a measure of local adaptation that accounts for limited dispersal
ℓₗᵢₘ = function (p)

    # m = max distance
    # s = step size (ie., resolution)
    # p = model parameters

    # unpack some parameters
    @unpack Bₕ,Bₚ = p

    # compute cross-covariance function
    C̄HP = C̄ₕₚ(p)
    CHP0 = Cₕₚ₀(p)

    ℓₕ = Bₕ * (C̄HP - CHP0)
    ℓₚ = Bₚ * (CHP0 - C̄HP)
    
    ℓ = [ℓₕ, ℓₚ]

    return(ℓ)

end

# returns a measure of local adaptation pretending dispersal follows an island model when it's actually limited
ℓₖₗₛ = function (p)

    # p = model parameters

    # unpack strengths of coevolutionary selection
    @unpack Bₕ,Bₚ = p
    
    ℓₕ = -Bₕ .* Cₕₚ₀(p)
    ℓₚ =  Bₚ .* Cₕₚ₀(p)

    ℓ = [ℓₕ, ℓₚ]

    return(ℓ)

end

# returns spatial covariance of trait means in the case of an island model
CₕₚISL = function (p)

    # p = model parameters

    # unpack strengths of coevolutionary selection
    @unpack Bₕ,Bₚ,Aₕ,Aₚ,Gₚ,Gₕ,ρₚ,ρₕ,σₚ,σₕ = p
    
    # use discretization of laplacian to justify analogous dispersal par σ^2

    num = Gₕ*Gₚ*( Bₚ*ρₚ*(σₚ+Gₚ*(Aₚ+Bₚ)) - Bₕ*ρₕ*(σₕ+Gₕ*(Aₕ-Bₕ)) )

    den₁ = 2*ρₕ*ρₚ*(σₕ+σₚ+Gₕ*(Aₕ-Bₕ)+Gₚ*(Aₚ+Bₚ))

    den₂ = σₕ*σₚ + σₕ*Gₚ*(Aₚ+Bₚ) + σₚ*Gₚ*(Aₕ-Bₕ) + Gₕ*Gₚ*(Aₕ*Aₚ+Aₕ*Bₚ-Aₚ*Bₕ)

    Cₕₚ = num/(den₁*den₂)

    return(Cₕₚ)

end

# returns a measure of local adaptation in the absence of gene-flow
ℓᵢₛₗ = function (p)

    # p = model parameters

    # unpack strengths of coevolutionary selection
    @unpack Bₕ,Bₚ = p
    
    ℓₕ = -Bₕ * CₕₚISL(p)
    ℓₚ =  Bₚ * CₕₚISL(p)

    ℓ = [ℓₕ, ℓₚ]

    return(ℓ)

end

# returns spatial covariance of trait means in the absence of gene-flow
CₕₚNGF = function (p)

    # p = model parameters

    # unpack strengths of coevolutionary selection
    @unpack Bₕ,Bₚ,Aₕ,Aₚ,Gₚ,Gₕ,ρₚ,ρₕ = p
    
    Cₕₚ = (Bₚ*(Aₚ+Bₚ)*Gₚ*ρₚ-Bₕ*(Aₕ-Bₕ)*Gₕ*ρₕ) / 
        ( 2*ρₕ*ρₚ*(Aₕ*Aₚ+Aₕ*Bₚ-Aₚ*Bₕ)*(Gₕ*(Aₕ-Bₕ)+Gₚ*(Aₚ+Bₚ)) )

    return(Cₕₚ)

end

# returns a measure of local adaptation in the absence of gene-flow
ℓₘ₌₀ = function (p)

    # p = model parameters

    # unpack strengths of coevolutionary selection
    @unpack Bₕ,Bₚ = p
    
    ℓₕ = -Bₕ .* CₕₚNGF(p)
    ℓₚ =  Bₚ .* CₕₚNGF(p)

    ℓ = [ℓₕ, ℓₚ]

    return(ℓ)

end

# index of coevolutionary advantage
𝓪 = function(p)

    @unpack Aₕ, Aₚ, Bₕ, Bₚ, Gₕ, Gₚ, ρₕ, ρₚ, σₕ, σₚ = p

    Vₕ = 1/(ρₕ*σₕ^2*(Aₕ-Bₕ))
    Vₚ = 1/(ρₚ*σₚ^2*(Aₚ+Bₚ))

    Vₕ⁰ = 1/(ρₕ*σₕ^2*Aₕ)
    Vₚ⁰ = 1/(ρₚ*σₚ^2*Aₚ)

end