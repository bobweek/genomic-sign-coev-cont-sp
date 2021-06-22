using LinearAlgebra
using Gadfly
using Alert
using DataFrames
using CSV

include("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/numerical-cov.jl")


# NOTATION

# 𝓛, this font is \bscrL

# ℓₛ = local adaptation of spp 𝑠
# 𝓖ₛ = global adaptation of spp 𝑠
# 𝓜ₛ = spatial average of growth rate for spp 𝑠
# 𝓜ₛₙ = spatial average of growth rate for spp 𝑠 with no biotic selection
# 𝓐ₛ = coev advantage of spp 𝑠

# 𝓜ₛ = ℓₛ + 𝓖ₛ
# 𝓐ₛ = 𝓜ₛ - 𝓜ₛₙ

#
# first part plots in (σₕ,σₚ)-coords
#

# background parameters
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, ρₕ = 10, ρₚ = 10, 𝓝ₕ=100, 𝓝ₚ=100, Aₕ = 1, Aₚ = 1, Bₕ = 0.01, Bₚ = 0.01, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)

# for plotting LA under the island model
σ = 10 .^ ((-1):0.01:2)
ISL_df₁₀ = DataFrame(σₕ = Float64[], σₚ = Float64[], LAₕ = Float64[], LAₚ = Float64[] )
for sₕ in σ
    for sₚ in σ

        p.σₕ = sₕ
        p.σₚ = sₚ
        # p.ρₕ = p.𝓝ₕ / (4*π*sₕ^2)
        # p.ρₚ = p.𝓝ₚ / (4*π*sₚ^2)
        LAᵢₛₗ = ℓᵢₛₗ(p) # LA under island model
        push!(ISL_df₁₀,(sₕ, sₚ, LAᵢₛₗ[1], LAᵢₛₗ[2]))

    end
end
ISLpl = Gadfly.plot(ISL_df₁₀,x=:σₕ, y=:σₚ, z=:LAₚ, Geom.contour(), Scale.x_log10, Scale.y_log10, Guide.xlabel("log₁₀(Host Dispersal Rate)"), Guide.ylabel("log₁₀(Parasite Dispersal Rate)"), Guide.Title("Parasite Local Adaptation Under Island Model"), Coord.Cartesian(ymax=1))
CSV.write("gits/genomic-sign-coev-cont-sp/phenotypic/julia/ISL.csv",ISL_df₁₀)

# for plotting LA under limited disp
σ = 10 .^ (0:0.01:3)
CLS_df₁₀₀ = DataFrame(σₕ = Float64[], σₚ = Float64[], LAₕ = Float64[], LAₚ = Float64[] )
LMD_df₁₀₀ = DataFrame(σₕ = Float64[], σₚ = Float64[], LAₕ = Float64[], LAₚ = Float64[] )
for sₕ in σ
    for sₚ in σ

        p.σₕ = sₕ
        p.σₚ = sₚ

        LAₖₗₛ = ℓₖₗₛ(p) # classical LA w limited disp
        LAₗᵢₘ = ℓₗᵢₘ(p) # modified LA w limited disp

        push!(CLS_df₁₀₀,(sₕ, sₚ, LAₖₗₛ[1], LAₖₗₛ[2]))
        push!(LMD_df₁₀₀,(sₕ, sₚ, LAₗᵢₘ[1], LAₗᵢₘ[2]))

    end
end
CLSpl = Gadfly.plot(CLS_df₁₀₀,x=:σₕ, y=:σₚ, z=:LAₚ, Geom.contour(), Scale.x_log10, Scale.y_log10, Guide.xlabel("log₁₀(Host Dispersal Distance)"), Guide.ylabel("log₁₀(Parasite Dispersal Distance)"), Guide.Title("Parasite Local Adaptation With Limited Dispersal"))
LMDpl = Gadfly.plot(LMD_df₁₀₀,x=:σₕ, y=:σₚ, z=:LAₚ, Geom.contour(), Scale.x_log10, Scale.y_log10, Guide.xlabel("log₁₀(Host Dispersal Distance)"), Guide.ylabel("log₁₀(Parasite Dispersal Distance)"), Guide.Title("Parasite Local Adaptation With Limited Dispersal"))

CSV.write("gits/genomic-sign-coev-cont-sp/phenotypic/julia/CLS.csv",CLS_df₁₀₀)
CSV.write("gits/genomic-sign-coev-cont-sp/phenotypic/julia/LMD.csv",LMD_df₁₀₀)

alert("Done")

#
# second part plots in (Bₕ,Bₚ)-coords
#

# resetting background parameters
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, ρₕ = 10, ρₚ = 10, Aₕ = 1, Aₚ = 1, Bₕ = 0.2, Bₚ = 0.2, σₕ = 5, σₚ = 5, rₕ = 0, rₚ = 0)


# parameters for numerical integration
m = √(π*(p.σₕ^2+p.σₚ^2)/2)
s = m/10

# for plotting LA under the island model
B = 0.001:0.002:0.1
ISL_df_B = DataFrame(Bₕ = Float64[], Bₚ = Float64[], LAₕ = Float64[], LAₚ = Float64[] )
for bₕ in B
    for bₚ in B

        p.Bₕ = bₕ
        p.Bₚ = bₚ
        LAᵢₛₗ = ℓᵢₛₗ(p) # LA under island model
        push!(ISL_df_B,(bₕ, bₚ, LAᵢₛₗ[1], LAᵢₛₗ[2]))

    end
end
ISLpl = Gadfly.plot(ISL_df_B,x=:Bₕ, y=:Bₚ, z=:LAₚ, Geom.contour())
CSV.write("gits/genomic-sign-coev-cont-sp/phenotypic/julia/ISL_B.csv",ISL_df_B)

# for plotting LA under limited disp
B = 0.001:0.002:0.1
CLS_df_B = DataFrame(Bₕ = Float64[], Bₚ = Float64[], LAₕ = Float64[], LAₚ = Float64[] )
LMD_df_B = DataFrame(Bₕ = Float64[], Bₚ = Float64[], LAₕ = Float64[], LAₚ = Float64[] )
for bₕ in B
    for bₚ in B

        p.Bₕ = bₕ
        p.Bₚ = bₚ

        LAₖₗₛ = ℓₖₗₛ(p) # classical LA w limited disp
        LAₗᵢₘ = ℓₗᵢₘ(m,s,p) # modified LA w limited disp

        push!(CLS_df_B,(bₕ, bₚ, LAₖₗₛ[1], LAₖₗₛ[2]))
        push!(LMD_df_B,(bₕ, bₚ, LAₗᵢₘ[1], LAₗᵢₘ[2]))

    end
end
Gadfly.plot(CLS_df_B, x=:Bₕ, y=:Bₚ, z=:LAₚ, Geom.contour(), Guide.Title("Parasite Local Adaptation With Limited Dispersal"))
Gadfly.plot(LMD_df_B, x=:Bₕ, y=:Bₚ, z=:LAₚ, Geom.contour(), Guide.Title("Parasite Local Adaptation With Limited Dispersal"))

CSV.write("gits/genomic-sign-coev-cont-sp/phenotypic/julia/CLS_B.csv",CLS_df_B)
CSV.write("gits/genomic-sign-coev-cont-sp/phenotypic/julia/LMD_B.csv",LMD_df_B)

alert("Done")