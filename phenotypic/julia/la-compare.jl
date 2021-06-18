using LinearAlgebra

include("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/numerical-cov.jl")

# range of disperal distances under consideration
σ = 0.1:1:10

CHP0 = zeros(length(σ),length(σ))
CHPd̄ = zeros(length(σ),length(σ))

iₕ = 1
for sₕ in σ
    iₚ = 1
    for sₚ in σ
        p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 1, Aₚ = 1, Bₕ = 0.2, Bₚ = 0.2, σₕ = sₕ, σₚ = sₚ, rₕ = 0, rₚ = 0)
        CHP0[iₕ,iₚ] = Cₕₚ₀(p)
        CHPd̄[iₕ,iₚ] = Cₕₚd̄(60,0.1,p)
        iₚ += 1
    end
    iₕ += 1
end

# range of coevolutionary selection under consideration
B = 0.001:0.002:0.01

CHP0 = zeros(length(B))
CHPd̄ = zeros(length(B))

i = 1
for b in B
    p = CoevPars(Gₚ = 10, Gₕ = 10, vₕ=10, vₚ=10, Nₕ = 10, Nₚ = 10, Aₕ = 0.1, Aₚ = 0.1, Bₕ = b, Bₚ = b, σₕ = 1, σₚ = 1, rₕ = 0, rₚ = 0)
    CHP0[i] = Cₕₚ₀(p)
    CHPd̄[i] = Cₕₚd̄(60,0.1,p)    
    i += 1
end

plot(B,CHP0)

plot(CHP0,CHPd̄)

plot(B,CHPd̄)
