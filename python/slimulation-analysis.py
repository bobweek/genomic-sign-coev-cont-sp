import pandas as pd
import numpy as np

ind0001 = pd.read_csv("python/ind-data/indData0011.csv")

hosts = ind0001.index[ind0001.spp==1]
Nh = len(hosts)
Np = len(ind0001.index[ind0001.spp==2])

traitpairs = []
paraperhost = []
hostedparas = []
for h in hosts:
    paraon = (ind0001.spp==2) & (ind0001.x==ind0001.x[h]) & (ind0001.y==ind0001.y[h])
    hostedparas.append(paraon)
    paraperhost.append(sum(paraon))
    for p in ind0001.index[paraon]:        
        traitpairs.append([ind0001["zh"][h], ind0001["zh"][p]]) # typo, zh should just be z, will fix later

# proportion of unhosted parasites
pr_unhosted = (Np - len(hostedparas))/Np

# proportion of unparasitized hosts
unparasitized = np.where(np.asarray(paraperhost)==0)[0]
pr_unparasitized = len(unparasitized)/Nh

# mean and variance of parasites per host
pprh_m = np.mean(paraperhost)
pprh_v = np.var(paraperhost)

np.corrcoef(np.transpose(traitpairs))

# host "selection matrix"
hselM = np.cov(np.transpose(ind0001[ind0001.spp==1][["zh","W"]]))
v_h = hselM[0,0] # phenotypic variance
beta_h = hselM[0,1]/v_h # global host selection gradient

# para "selection matrix"
pselM = np.cov(np.transpose(ind0001[ind0001.spp==2][["zh","W"]]))
v_p = pselM[0,0] # phenotypic variance
beta_p = pselM[0,1]/v_p # global parasite selection gradient

# discr sp by intxn dist and compute fraction occupied by each spp, mean and var of abund per square for each spp, mean and var of selection gradient per square for each spp

# then iterate thru each recorded tick to accumulate time-series for each variable of interest