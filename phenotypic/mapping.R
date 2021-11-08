host_df = read.csv("julia/testing-mapping/host_df-1-3.csv");
para_df = read.csv("julia/testing-mapping/para_df-1-1.csv");

host_df$trait

# to get rho_S, drop a few disks of radius R_S
# make sure they don't intersect
# then estimate rho_S by average density of inds in each disk

ndisk = 20
RH = 0.005	# competition radii
RP = 0.005

# pick their centers, calculate distance & count individuals
centersH = matrix(runif(2*ndisk, min=RH, max=1-RH),ndisk,2)
centersP = matrix(runif(2*ndisk, min=RP, max=1-RP),ndisk,2)
distsH = matrix(0,ndisk,ndisk)
distsP = matrix(0,ndisk,ndisk)
for(i in 1:ndisk){
  for(j in 1:ndisk){
    distsH[i,j] = sqrt( sum((centersH[i,]-centersH[j,])^2) )
    distsP[i,j] = sqrt( sum((centersP[i,]-centersP[j,])^2) )
  }
}
indsH = list()
indsP = list()
nindsH = rep(0,3)
nindsP = rep(0,3)
for(i in 1:ndisk){
  indsH[[i]] = which(rowSums((cbind(host_df1$x1,host_df1$x2)-centersH[i,])^2)<RH)
  indsP[[i]] = which(rowSums((cbind(para_df1$x1,para_df1$x2)-centersP[i,])^2)<RP)
  nindsH[i] = length(indsH[[i]])
  nindsP[i] = length(indsP[[i]])
}

# ensures they're not overlapping, not empty and not too close to the edge
mH = min(distsH[upper.tri(distsH)])
mP = min(distsP[upper.tri(distsP)])
while(mH<(2*RH) | min(nindsH)==0){
  centersH = matrix(runif(2*ndisk, min=RH, max=1-RH),ndisk,2)
  for(i in 1:ndisk){
    for(j in 1:ndisk){
      distsH[i,j] = sqrt( sum((centersH[i,]-centersH[j,])^2) )
    }
  }
  mH = min(distsH[upper.tri(distsH)])
  for(i in 1:ndisk){
    indsH[[i]] = which(rowSums((cbind(host_df1$x1,host_df1$x2)-centersH[i,])^2)<RH)
    nindsH[i] = length(indsH[[i]])
  }
}
while(mP<(2*RP) | min(nindsP)==0){
  centersP = matrix(runif(2*ndisk, min=RP, max=1-RP),ndisk,2)
  for(i in 1:ndisk){
    for(j in 1:ndisk){
      distsP[i,j] = sqrt( sum((centersP[i,]-centersP[j,])^2) )
    }
  }
  mP = min(distsP[upper.tri(distsP)])
  for(i in 1:ndisk){
    indsP[[i]] = which(rowSums((cbind(para_df1$x1,para_df1$x2)-centersP[i,])^2)<RP)
    nindsP[i] = length(indsP[[i]])
  }
}

zbarH_vec = rep(0,ndisk)
zbarP_vec = rep(0,ndisk)
GH_vec = rep(0,ndisk)
GP_vec = rep(0,ndisk)
for(i in 1:ndisk){
  zbarH_vec[i] = mean(host_df1$trait[indsH[[i]]])
  zbarP_vec[i] = mean(para_df1$trait[indsP[[i]]])
  nH = length(indsH[[i]])
  nP = length(indsP[[i]])
  GH_vec[i] = (nH-1)*var(host_df1$trait[indsH[[i]]])/nH
  GP_vec[i] = (nP-1)*var(para_df1$trait[indsP[[i]]])/nP
}

GH_hat = mean(GH_vec)
GP_hat = mean(GP_vec)

pi = 3.14159

rhoH_hat = mean(nindsH) #/(pi*RH^2)
rhoP_hat = mean(nindsP) #/(pi*RP^2)

rhoH_var = var(nindsH) #/(pi*RH^2)
rhoP_var = var(nindsP) #/(pi*RP^2)

gamma = 0.1
pimax = 1.0
Ri = 0.01

iotaP = 1.1
alphP = 1.2
kappaP = 0.98
muP = 0.1
EP = 0.001
AP = 0.002

iotaH = 0.99
alphH = 1.4
kappaH = 0.975
EH = 0.001
muH = 0.1
AH = 0.002

minthis = function(X){
  
  rhoH = X[1]
  rhoP = X[2]
  BP = gamma*pimax*(iotaP-1)*(1-exp(-2*pi*Ri^2*rhoH))
  vP = EP + sqrt(muP/(AP+BP))
  BH = -gamma*pimax*(iotaH-1)*rhoP/rhoH 
  
  rP = log(alphP)+pimax*(iotaP-1)*(1-exp(-2*pi*Ri^2*rhoH))
  rH = log(alphH)+pimax*(iotaH-1)*(1-gamma*vP/2)*rhoP/rhoH
  
  out = c((rhoH+(rH-0.5*sqrt(muH*abs(AH-BH)))/log(kappaH))^2,
          (rhoP+(rP-0.5*sqrt(muP*(AP+BP)))/log(kappaP))^2 )
  
  return(sqrt(sum(out)))
}

num_rho = optim(rexp(2,1/10),minthis)

rhoH = num_rho$par[1] #/(2*pi*RH^2)
rhoP = num_rho$par[2] #/(2*pi*RP^2)

BP = gamma*pimax*(iotaP-1)*(1-exp(-2*pi*Ri^2*rhoH))
vP = EP + sqrt(muP/(AP+BP))
BH = -gamma*pimax*(iotaH-1)*rhoP/rhoH 
vH = EH + sqrt(muH/(AH-BH))
rP = log(alphP)+pimax*(iotaP-1)*(1-exp(-2*pi*Ri^2*rhoH))
rH = log(alphH)+pimax*(iotaH-1)*(1-gamma*vP/2)*rhoP/rhoH
GH = sqrt(muH/(AH-BH))
GP = sqrt(muP/(AP+BP))


ExpObs.df = data.frame(Parameter=c("host density","parasite density", "host G", "parasite G"), Expectation=c(rhoH, rhoP, GH, GP), Observation=c(rhoH_hat, rhoP_hat, GH_hat, GP_hat))