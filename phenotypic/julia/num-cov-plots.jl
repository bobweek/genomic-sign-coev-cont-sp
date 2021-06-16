include("/home/bb/gits/white.noise.community.ecology/ibm_functions_structs.jl")

# same dispersal distances
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
eqLA = plotLocAdapt(60,0.1,p,"σ")
title!("")
# plot(eqCorr, eqLA, xlims=(0,40), layout = (2,1), size = (400,600))
plot(eqLA, xlims=(0,40), title = "", size = (500,400))
savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/eq-pars.png")

# parasite disperses further than host
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, σₕ = 10, σₚ = 100, rₕ = 0, rₚ = 0)
# Corrp = plotSpCorr(60,0.1,p)
LAp = plotLocAdapt(60,0.1,p,"σ")
# plot(Corrp, LAp, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/p-further.png")


# host disperses a tiny bit further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, σₕ = 10.5, σₚ = 10, rₕ = 0, rₚ = 0)
# Corrhl = plotSpCorr(60,0.1,p)
LAht = plotLocAdapt(60,0.1,p,"σ")
# plot(Corrhl, LAhl, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-l-further.png")

# host disperses a little bit further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, σₕ = 20, σₚ = 10, rₕ = 0, rₚ = 0)
# Corrhl = plotSpCorr(60,0.1,p)
LAhl = plotLocAdapt(60,0.1,p,"σ")
# plot(Corrhl, LAhl, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-l-further.png")

# host disperses an intermediate bit further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, σₕ = 50, σₚ = 10, rₕ = 0, rₚ = 0)
# Corrhi = plotSpCorr(60,0.1,p)
LAhi = plotLocAdapt(60,0.1,p,"σ")
# plot(Corrhi, LAhi, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-i-further.png")

# host disperses much further than parasite
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, σₕ = 100, σₚ = 10, rₕ = 0, rₚ = 0)
# Corrh = plotSpCorr(60,0.1,p)
LAh = plotLocAdapt(60,0.1,p,"σ")
# plot(Corrh, LAh, xlims=(0,40), layout = (2,1), size = (400,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/h-further.png")

# plot(Corrp, eqCorr, Corrhl, Corrhi, Corrh, eqLA, LAp, LAhl, LAhi, LAh, xlims=(0,40), layout=(2,5), size=(1600,600))
# savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/corrs-las.png")

# equal strength
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.01, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
eqBLA = plotLocAdapt(60,0.1,p,"B")

# host stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.02, Bₚ = 0.01, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
hsLA = plotLocAdapt(60,0.1,p,"B")

# host little stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.0106, Bₚ = 0.01, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
hlsLA = plotLocAdapt(60,0.1,p,"B")

# host much stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.1, Bₚ = 0.01, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
hmsLA = plotLocAdapt(60,0.1,p,"B")

# parasite stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.02, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
psLA = plotLocAdapt(60,0.1,p,"B")

# parasite much stronger
p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.2, Aₚ = 0.2, Bₕ = 0.01, Bₚ = 0.1, σₕ = 10, σₚ = 10, rₕ = 0, rₚ = 0)
# eqCorr = plotSpCorr(60,0.1,p)
pmsLA = plotLocAdapt(60,0.1,p,"B")

plot(LAp, eqLA, LAht, LAhl, LAhi, LAh, xlims=(0,40), size=(1200,600))
savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/d-las.png")

plot(pmsLA, eqBLA, hlsLA, hmsLA, xlims=(0,40), size=(1200,600))
savefig("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/b-las.png")
