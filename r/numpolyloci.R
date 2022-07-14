# 2 Ne lineages
# mu is pr mutation per step
# th = 4 Ne mu is population-scaled mutation rate
# but mutations actually occur at rate th/2 per lineage
# so whole population mutation rate is 2*Ne*th/2=4*mu*Ne^2

Ne = 1e2
mu = 1/Ne

reps = c()
for(i in 1:100){
  m = cumsum(rexp(20*Ne,4*mu*Ne^2)) # mutation times
  ply = c()
  for(l in m){
    t = c()
    for(i in (2*Ne):2){ # draw coal times
      t = c(t,rexp(1,choose(i,2)))
    }
    t = cumsum(t)
    wch = which(t<l)
    if(length(wch)>0){
      p = (2*Ne-max(wch)-1+2)/(2*Ne)
    } else p = 1
    ply = c(ply,rbinom(1,1,p))
  }
  reps = c(reps,sum(ply))
}
mean(reps)

g = 0.577215664
Sw = 4*Ne*mu*(g+log(2*Ne))
