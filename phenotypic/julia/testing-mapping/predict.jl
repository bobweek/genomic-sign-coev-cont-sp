#
# trying to guess simulated results...
#

# be sure to cd into the folder containing this file first
include("../ibm_functions_structs.jl")

# values of competition effects to check
κₕ = 0.9;
κₚ = 0.9;

# mutation rates to check
μₕ = 0.1;
μₚ = 0.1;

# solving for expected densities
function ρsolve(X)
    ρ̃ₕ = exp(X[1])
    ρ̃ₚ = exp(X[2])


    B̃ₕ = Bₕ(γ, πₘ, ιₕ, ρ̃ₚ, ρ̃ₕ)
    B̃ₚ = Bₚ(γ, πₘ, ιₚ, Nᵨ(ρ̃ₕ,Rᵢ))
    ṽₕ = v(Eₕ, Gₚ(μₕ, Aₕ, B̃ₕ))
    ṽₚ = v(Eₚ, Gₚ(μₚ, Aₚ, B̃ₚ))
    r̃ₕ = rₕ(αₕ, Aₕ, ṽₕ, πₘ, ιₕ, γ, ṽₚ, ρ̃ₚ, ρ̃ₕ, B̃ₕ)
    r̃ₚ = rₚ(αₚ, Aₚ, ṽₚ, πₘ, ιₚ, Nᵨ(ρ̃ₕ,Rᵢ), B̃ₚ, ṽₕ)
    z̃ₕ = z̄ₕ(Aₕ, θ₀ₕ, B̃ₕ, Aₚ, θ₀ₚ, B̃ₚ)
    z̃ₚ = z̄ₚ(Aₚ, θ₀ₚ, B̃ₚ, Aₕ, θ₀ₕ, B̃ₕ)
    m̃ₕ = mₕ(r̃ₕ, z̃ₕ, Aₕ, θ₀ₕ, B̃ₕ, z̃ₚ)
    m̃ₚ = mₚ(r̃ₚ, z̃ₚ, Aₚ, θ₀ₚ, B̃ₚ, z̃ₕ)
    Ñₕ = Nₕ(κₕ, m̃ₕ, μₕ, Aₕ, B̃ₕ)
    Ñₚ = Nₚ(κₚ, m̃ₚ, μₚ, Aₚ, B̃ₚ)

    val = √((ρ̃ₕ - ρ(Ñₕ,Rₕ))^2 + (ρ̃ₚ - ρ(Ñₚ,Rₚ))^2)

    return val

end

# for checking exp densities, coll vars and char lngths
ρ̃ₕ = 0.0
ρ̃ₚ = 0.0
B̃ₕ = 0.0
B̃ₚ = 0.0
ṽₕ = 0.0
ṽₚ = 0.0
r̃ₕ = 0.0
r̃ₚ = 0.0
z̃ₕ = 0.0
z̃ₚ = 0.0
m̃ₕ = 0.0
m̃ₚ = 0.0
Ñₕ = 0.0
Ñₚ = 0.0
Ṽₕ = 0.0
Ṽₚ = 0.0
ξ̃ₕ = 0.0
ξ̃ₚ = 0.0
keep = false
while !keep
    initp = rand(Uniform(1, 20), 2)    
    try
        res = optimize(ρsolve, initp, Optim.Options(iterations = Int64(1e5)))
        ρ̃ₕ = exp(Optim.minimizer(res)[1])
        ρ̃ₚ = exp(Optim.minimizer(res)[2])
        keep = true
    catch err
        if isa(err, DomainError)
            print("\n")
            print("poopy!")
            print("\n")
        end
    end
end
ρ̃ₕ
ρ̃ₚ
