host_df1 = read.csv("julia/ibm/host_df1.csv");
para_df1 = read.csv("julia/ibm/para_df1.csv");


host_N = length(host_df1$id)
para_N = length(para_df1$id)
host_gdist = matrix(0,host_N,host_N)
host_pdist = matrix(0,host_N,host_N)
for(i in 1:host_N){
  for(j in 1:host_N){
    xi = c(host_df1$x1[i], host_df1$x2[i])
    xj = c(host_df1$x1[j], host_df1$x2[j])
    zi = host_df1$trait[i]
    zj = host_df1$trait[j]
    host_gdist[i,j] = sqrt( sum((xi-xj)^2) )
    host_pdist[i,j] = abs(zi-zj)
  }
}
para_gdist = matrix(0,para_N,para_N)
para_pdist = matrix(0,para_N,para_N)
for(i in 1:para_N){
  for(j in 1:para_N){
    xi = c(para_df1$x1[i], para_df1$x2[i])
    xj = c(para_df1$x1[j], para_df1$x2[j])
    zi = para_df1$trait[i]
    zj = para_df1$trait[j]
    para_gdist[i,j] = sqrt( sum((xi-xj)^2) )
    para_pdist[i,j] = abs(zi-zj)
  }
}
hp_gdist = matrix(0,host_N,para_N)
hp_pdist = matrix(0,host_N,para_N)
for(i in 1:host_N){
  for(j in 1:para_N){
    xi = c(host_df1$x1[i], host_df1$x2[i])
    xj = c(para_df1$x1[j], para_df1$x2[j])
    zi = host_df1$trait[i]
    zj = para_df1$trait[j]
    hp_gdist[i,j] = sqrt( sum((xi-xj)^2) )
    hp_pdist[i,j] = abs(zi-zj)
  }
}

par(mfrow=c(1,3))

sample_ids = sample(1:choose(host_N,2),1000)
hgeo = host_gdist[upper.tri(host_gdist)][sample_ids]
hpheno = host_pdist[upper.tri(host_pdist)][sample_ids]
plot(hgeo,hpheno,xlab="Geographic Distance",ylab="Phenotypic Difference")
text(1,max(hpheno)-1,"Host")

sample_ids = sample(1:choose(para_N,2),1000)
pgeo = para_gdist[upper.tri(para_gdist)][sample_ids]
ppheno = para_pdist[upper.tri(para_pdist)][sample_ids]
plot(pgeo,ppheno,xlab="Geographic Distance",ylab="Phenotypic Difference")
text(1,max(ppheno)-1,"Parasite")

hpgeo = hp_gdist[upper.tri(hp_gdist)][sample_ids]
hppheno = hp_pdist[upper.tri(hp_pdist)][sample_ids]
plot(hpgeo,hppheno,xlab="Geographic Distance",ylab="Phenotypic Difference")
text(1,max(ppheno)-1,"Host-Parasite")
