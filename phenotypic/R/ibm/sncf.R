#
# this script is for analyzing spatial correlations of simulated trait data
#

# this is the package used to get non-parametric corr fct
require(ncf)

# load in host/parasite trait/space data
host_df1 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/host_df1.csv");
para_df1 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/para_df1.csv");
host_df2 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/host_df2.csv");
para_df2 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/para_df2.csv");
host_df3 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/host_df3.csv");
para_df3 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/para_df3.csv");
host_df4 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/host_df4.csv");
para_df4 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/para_df4.csv");
host_df5 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/host_df5.csv");
para_df5 = read.csv("/home/bb/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ibm/para_df5.csv");

host_df = list(host_df1,host_df2,host_df3,host_df4,host_df5)
para_df = list(para_df1,para_df2,para_df3,para_df4,para_df5)

# discretize space into S^2 units
S = 20;
xbins = (1:S)/S;

hmeans = list(NULL,NULL,NULL,NULL,NULL) # format = c(x1 bin, x2 bin, trait mean)
pmeans = list(NULL,NULL,NULL,NULL,NULL)
tmeans = list(NULL,NULL,NULL,NULL,NULL) # format = c(x1 bin, x2 bin, h trait mean, p trait mean)

for(k in 1:5){
  # denotes which horiz/vert bin each ind belongs to
  hbins = matrix(0,nrow=length(host_df[[k]]$id),ncol=2);
  pbins = matrix(0,nrow=length(para_df[[k]]$id),ncol=2);
  
  # first bins
  hsbst = subset(host_df[[k]], x1<xbins[1])
  psbst = subset(para_df[[k]], x1<xbins[1])
  hbins[hsbst$id,1] = 1
  pbins[psbst$id,1] = 1
  
  hsbst = subset(host_df[[k]], x2<xbins[1])
  psbst = subset(para_df[[k]], x2<xbins[1])
  hbins[hsbst$id,2] = 1
  pbins[psbst$id,2] = 1
  
  # rest of the bins
  for(i in 2:S){
    
    hsbst = subset(host_df[[k]], x1>xbins[i-1] & x1<xbins[i])
    psbst = subset(para_df[[k]], x1>xbins[i-1] & x1<xbins[i])
    hbins[hsbst$id,1] = i
    pbins[psbst$id,1] = i
    
    hsbst = subset(host_df[[k]], x2>xbins[i-1] & x2<xbins[i])
    psbst = subset(para_df[[k]], x2>xbins[i-1] & x2<xbins[i])
    hbins[hsbst$id,2] = i
    pbins[psbst$id,2] = i
    
  }
  
  # accumulate average trait values for each non-empty bin
  for(i in 1:S){
    for(j in 1:S){
      
      hi = which(hbins[,1]==i)
      hj = which(hbins[,2]==j)
      hinds = intersect(hi,hj)
      if(length(hinds)>0) hmeans[[k]] = rbind(hmeans[[k]],c(i/S,j/S,mean(host_df[[k]]$trait[hinds])))
      
      pi = which(pbins[,1]==i)
      pj = which(pbins[,2]==j)
      pinds = intersect(pi,pj)
      if(length(pinds)>0) pmeans[[k]] = rbind(pmeans[[k]],c(i/S,j/S,mean(para_df[[k]]$trait[pinds])))
      
      if(length(hinds)>0 & length(pinds)>0) tmeans[[k]] = rbind(tmeans[[k]],c(i/S,j/S,mean(host_df[[k]]$trait[hinds]),mean(para_df[[k]]$trait[pinds])))
      
    }
  }
}

# now we need to go through and lump together replicates at each location
# i think we'll need same number of replicates for each location
intsct = list(c(),c(),c(),c(),c())
for(i in 1:S){
  for(j in 1:S){
    
    # rw contains row numbers corresponding to (i,j)th bin
    rw = c(0,0,0,0,0)
    for(k in 1:5){
      if(length(which(tmeans[[k]][,1]==i/S & tmeans[[k]][,2]==j/S))>0){
        rw[k] = which(tmeans[[k]][,1]==i/S & tmeans[[k]][,2]==j/S)
      }
    }

    if( min(rw)>0 ){
      print('poopi')
      for(k in 1:5) intsct[[k]] = rbind(intsct[[k]], tmeans[[k]][rw[k],])
    }

  }
}



fit <- Sncf.srf(x = tmeans[,1], y = tmeans[,2], z = tmeans[,3:4], avg = NULL, resamp = 0) 

fit <- Sncf.srf(x = hmeans[,1], y = hmeans[,2], z = hmeans[,3], avg = NULL, corr = TRUE, resamp = 0) 

fit <- Sncf.srf(x = tmeans[,1], y = tmeans[,2], z = tmeans[,3], w = tmeans[,4], avg = NULL, avg2 = NULL, corr = TRUE, resamp = 0) 

fit <- Sncf(x = tmeans[,1], y = tmeans[,2], z = tmeans[,3:4], resamp = 0)

cor(tmeans[[4]][,3:4])

plot(fit)
summary(fit)
