using Optim, Distances, Distributions, LinearAlgebra, SpecialFunctions

# making the matern in 
# my preferred parameterization
function C(d, V, ξ)
    if (d/ξ) > 1e-300 && (d/ξ) < 1e+5
        try
            val = V * (√2 * d / ξ) * besselk(1.0, √2 * d / ξ)                            
        catch err
            # if evaluating the besselk function returns an error,
            # this will spit out the variables used to generate the error
            if isa(err, SpecialFunctions.AmosException)
                print("\n")
                print([d,V,ξ])
                print("\n")                
            end
        end
        
    elseif (d/ξ) > 1e+5
        val = 0
    else
        val = V
    end
    return val
end

# spits out covariance matrix given locations
function Cmat(x, V, ξ, v)
    n = size(x)[2]
    val = zeros(n, n)
    for i = 1:n
        for j = 1:n
            dd = euclidean(x[:, i], x[:, j])
            val[i, j] = C(dd, V, ξ)
        end
    end
    return val + v * Matrix(I, n, n)
end

# provides MLE of expectation of local mean trait
function hattildez(z, S)
    val = 0.0
    n = length(z)
    Sinv = inv(S)
    for i = 1:n
        for j = 1:n
            val = val + Sinv[i, j] * (z[i] + z[j]) / 2
        end
    end
    return (val / sum(Sinv))
end

# log-likelihood function
function ell(z, x, V, ξ, v)
    n = length(z)
    S = Cmat(x, V, ξ, v)
    ztilde = hattildez(z, S)
    return log(det(S)) + transpose(z .- ztilde) * inv(S) * (z .- ztilde)
end
