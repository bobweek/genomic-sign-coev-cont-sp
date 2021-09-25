include("mle-ftcs.jl")

# sample size of simulations for each par combo
N = 10

# number of individuals measured
n = [50 100]

# char scale
ξ = [0.01 0.1 1]

# marginal var
V = [0.1 1 10]

# within pop var
v = [0.1 1 10]

ncombo = length(v) * length(V) * length(ξ) * length(n)

# mean
ztilde = 0

# containers for input and estim pars
input = fill(0.0, 3, ncombo, N)
mle = fill(0.0, 3, ncombo, N)

for j = 1:N
    i = 1
    for nn in n

        # their locations, uniform on unit square
        X = rand(Uniform(), 2, nn)

        for ξξ in ξ
            for VV in V
                for vv in v

                    # covariance matrix
                    S = Cmat(X, VV, ξξ, vv)

                    cnvrg = false
                    while !cnvrg

                        # trait values
                        Z = rand(MvNormal(fill(ztilde, nn), S))

                        # function to minimize  
                        function minThis(Y)
                            VVV = exp(Y[1])
                            ξξξ = exp(Y[2])
                            vvv = exp(Y[3])
                            val = ell(Z, X, VVV, ξξξ, vvv)
                            return val
                        end

                        # do the optimization
                        res = optimize(minThis, rand(Exponential(1), 3))

                        # if the optimizer converged
                        # save the results
                        cnvrg = Optim.converged(res)
                        if cnvrg
                            mle[:, i, j] = exp.(Optim.minimizer(res))
                            input[:, i, j] = [VV, ξξ, vv]
                            i += 1
                        end
                    end
                end
            end
        end
    end
end

# first 27 of mle corresponds to nn=50, second 27 to nn=1000

# exporting results to csv
using DataFrames, CSV

ids = []
Vs = []
ξs = []
vs = []
for j = 1:ncombo
    for i = 1:N
      append!(ids, j)      
    end
    append!(Vs, mle[1,j,:])
    append!(ξs, mle[2,j,:])
    append!(vs, mle[3,j,:])
end

# first column labels the parameter combo, this id is repeated N times for each combo
# second column is estimates of V
# third is estimates of ξ
# fourth is v
mle_df = DataFrame(id = ids, V = Vs, ξ = ξs, v = vs);

# dataframe for input parameters
input_df = DataFrame(id = collect(1:ncombo), V = input[1,:,1], ξ = input[2,:,1], v = input[3,:,1] );

CSV.write("mle.csv",mle_df)
CSV.write("input.csv",input_df)
