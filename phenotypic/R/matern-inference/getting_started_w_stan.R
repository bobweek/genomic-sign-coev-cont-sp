################################################################
################################################################
#	script for getting started in STAN
################################################################
################################################################

# This is the stan model block
#
# the reference manual can be found at 
# http://mc-stan.org/users/documentation/
#
# basically though, it's got 6 components:
# 1) functions block where you can define your own useful functions
# 2) data block that specifies what data you're passing to the model
# 3) transformed data block for any data transformations you want to 
#	 make in the model (only evaluated once)
# 4) parameters block where you specify what parameters you're inferring
# 5) transformed parameters block that allows you to specify parameter 
#	 transforms (re-evaluated when parameter values are updated during 
#	 inference)
# 6) model block where you specify priors and likelihood

stanBlock <- 
"functions {
	matrix spCov(int N, real V, real xi, matrix D) {
		matrix[N,N] modDist;
		matrix[N,N] parCov;
		modDist = sqrt(2) * D / xi;
		for(i in 1:N){
			for(j in i:N){
				parCov[i,j] =  D[i,j] == 0 ? V : V * modDist[i,j] * modified_bessel_second_kind(1,modDist[i,j]);
				parCov[j,i] = parCov[i,j];
			}
		}
		return parCov;	
	}
}
data {
	int<lower=2> N; 	  			// number of samples
	matrix[N, N] geoDist; 			// matrix of pairwise geographic distance 
	vector[N] P;					// phenotype measured in each individual
}
parameters {
	real<lower=0> V;				// collocated variance
	real<lower=0> xi;				// characteristic length
	real zbar;						// mean phenotype
}
transformed parameters {
	matrix[N,N] parCov;				// this specifies the parametric covariance matrix
	vector[N] zbarVec;					// vectorized phenotype mean
	parCov = spCov(N, V, xi, geoDist);
	zbarVec = rep_vector(zbar,N);
}
model {
	V ~ exp(1);										// prior on V
	xi ~ exp(1);										// prior on xi
	zbar ~ normal(0,1);										// prior on zbar
	P ~ zbarlti_normal(zbarVec,parCov);							// prior on global covariance
}
"

library(rstan)

# compile the model here
#
# note that every time you update the model code
#	you'll have to recompile, which takes some time
#
# also note that this might spit a bunch of warnings,
#	which can usually be ignored
myMod <- stan_model(model_code=stanBlock)

################################
#	sizbarlate a spatial process
#		for N populations
################################

N <- 100 #number of individuals
coords <- cbind(runif(N),runif(N))
geoDist <- fields::rdist(coords)
params <- list("V" = 1.3,
				"xi" = 0.7,
				"zbar" = 1.1)
simCov <- sqrt(2)*params$V*(geoDist/params$xi)*besselK(x=sqrt(2)*(geoDist/params$xi),nu=1)
diag(simCov) <- params$V
simPheno <- MASS::mvrnorm(n=1,zbar=rep(params$zbar,N),Sigma=simCov)


################################
#	use the stan model to infer 
#		sizbarlating parameters
################################

dataBlock <- list("N" = N,
				 "geoDist" = geoDist,
				 "P" = simPheno)
fit <- sampling(object = myMod,
				 data = dataBlock,iter=2e3,
				 chains = 2)

# check and see how it looks
V <- extract(fit,"V",perzbarte=FALSE,inc_warzbarp=FALSE)[,,1]
xi <- extract(fit,"xi",perzbarte=FALSE,inc_warzbarp=FALSE)[,,1]
zbar <- extract(fit,"zbar",perzbarte=FALSE,inc_warzbarp=FALSE)[,,1]

par(mfrow=c(1,3))
matplot(V,type='l',col=adjustcolor("black",0.5),ylim=range(c(V,params$V)))
	abline(h=params$V,col="red",lty=2)
matplot(xi,type='l',col=adjustcolor("black",0.5),ylim=range(c(xi,params$xi)))
	abline(h=params$xi,col="red",lty=2)
matplot(zbar,type='l',col=adjustcolor("black",0.5),ylim=range(c(zbar,params$zbar)))
	abline(h=params$zbar,col="red",lty=2)
