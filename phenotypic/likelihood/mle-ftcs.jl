using Optim, Distances, Distributions, SpecialFunctions, LinearAlgebra

# for coloring data
function logit(x)
    return 1 / (1 + exp(-x))
end

# making the matern in 
# my preferred parameterization
function C(d, V, ξ)
    if d > 0
        val = V * (√2 * d / ξ) * besselk(1, √2 * d / ξ)
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
    val = 0
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
