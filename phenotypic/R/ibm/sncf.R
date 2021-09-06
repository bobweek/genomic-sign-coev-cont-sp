#
# this script is for analyzing spatial correlations of simulated trait data
#

# this is the package used to get non-parametric corr fct
require(ncf)

# load in host/parasite trait/space data
host_df = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/host_df.csv");
para_df = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/para_df.csv");

# discretize space into S^2 units
S = 20;
xbins = (1:S)/S;

# denotes which horiz/vert bin each ind belongs to
hbins = matrix(0,nrow=length(host_df$id),ncol=2);
pbins = matrix(0,nrow=length(para_df$id),ncol=2);

# first bins
hsbst = subset(host_df, x1<xbins[1])
psbst = subset(para_df, x1<xbins[1])
hbins[hsbst$id,1] = 1
pbins[psbst$id,1] = 1

hsbst = subset(host_df, x2<xbins[1])
psbst = subset(para_df, x2<xbins[1])
hbins[hsbst$id,2] = 1
pbins[psbst$id,2] = 1

# rest of the bins
for(i in 2:S){
  
  hsbst = subset(host_df, x1>xbins[i-1] & x1<xbins[i])
  psbst = subset(para_df, x1>xbins[i-1] & x1<xbins[i])
  hbins[hsbst$id,1] = i
  pbins[psbst$id,1] = i
  
  hsbst = subset(host_df, x2>xbins[i-1] & x2<xbins[i])
  psbst = subset(para_df, x2>xbins[i-1] & x2<xbins[i])
  hbins[hsbst$id,2] = i
  pbins[psbst$id,2] = i
  
}

# accumulate average trait values for each non-empty bin
hmeans = NULL # format = c(x1 bin, x2 bin, trait mean)
pmeans = NULL
tmeans = NULL # format = c(x1 bin, x2 bin, h trait mean, p trait mean)
for(i in 1:S){
  for(j in 1:S){
    
    hi = which(hbins[,1]==i)
    hj = which(hbins[,2]==j)
    hinds = intersect(hi,hj)
    if(length(hinds)>0) hmeans = rbind(hmeans,c(i/S,j/S,mean(host_df$trait[hinds])))
    
    pi = which(pbins[,1]==i)
    pj = which(pbins[,2]==j)
    pinds = intersect(pi,pj)
    if(length(pinds)>0) pmeans = rbind(pmeans,c(i/S,j/S,mean(para_df$trait[pinds])))
    
    if(length(hinds)>0 & length(pinds)>0) tmeans = rbind(tmeans,c(i/S,j/S,mean(host_df$trait[hinds]),mean(para_df$trait[pinds])))
    
  }
}

fit <- Sncf.srf(x = tmeans[,1], y = tmeans[,2], z = tmeans[,3:4], avg = NULL, resamp = 0) 

fit <- Sncf.srf(x = hmeans[,1], y = hmeans[,2], z = hmeans[,3], avg = NULL, corr = TRUE, resamp = 0) 

fit <- Sncf.srf(x = tmeans[,1], y = tmeans[,2], z = tmeans[,3], w = tmeans[,4], avg = NULL, avg2 = NULL, corr = TRUE, resamp = 0) 

cor(tmeans[,3:4])

fit <- Sncf(x = tmeans[,1], y = tmeans[,2], z = tmeans[,3:4], resamp = 0)

plot(fit)
summary(fit)
