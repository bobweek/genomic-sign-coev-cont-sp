import pandas as pd
import numpy as np
import os
import glob

csvCounter = len(glob.glob(os.path.expanduser("~/gsccs-data/ind-data/*.csv")))-1

# selection parameters
pfname = "~/gsccs-data/params.csv"
pars = pd.read_csv(pfname)
γ  = pars["γ"]
sₕ = pars["sₕ"]
sₚ = pars["sₚ"]

summary_time_series = pd.DataFrame(columns=[
            "Nh","Np","zh","zp","vh","vp","pr_unhosted","pr_unparasitized","pprh_m","pprh_v","pntcorr"])

txt = "{time:d}"
time_pts = np.arange(0,csvCounter)
for t in time_pts:

    fname = "~/gsccs-data/ind-data/indData"+txt.format(time = t).zfill(4)+".csv"
    inds = pd.read_csv(fname)

    hosts = inds.index[inds.spp==1]
    paras = inds.index[inds.spp==2]

    zh = inds.z[inds.spp==1]
    zp = inds.z[inds.spp==2]

    # Bₕ = np.zeros(len(inds.x))
    # Bₚ = np.zeros(len(inds.x)) + 1 # any unhosted parasite has unity biotic component of fitness
    traitpairs = []
    paraperhost = []
    hostedparas = np.zeros(len(inds.x))
    for h in hosts:
        paraon = (inds.spp==2) & (inds.x==inds.x[h]) & (inds.y==inds.y[h])
        paraperhost.append(sum(paraon))
        # BH = 1
        for p in inds.index[paraon]:
            hostedparas[p] = 1
            traitpairs.append([zh[h],zp[p]])
        #     a = np.exp(-γ*abs(zₕ[h]-zₚ[p]))
        #     Bₚ[p] = np.exp(sₚ)*a + (1-a)
        #     BH *= np.exp(sₕ)*a + (1-a)
        # Bₕ[h] = BH

    # global census sizes
    Nh = len(hosts)    
    Np = len(paras)

    # global trait variances
    vh = np.var(zh)
    vp = np.var(zp)

    # global mean traits
    zh = np.mean(zh)
    zp = np.mean(zp)

    # proportion of unhosted parasites
    pr_unhosted = (Np - sum(hostedparas))/Np

    # proportion of unparasitized hosts
    unparasitized = np.where(np.asarray(paraperhost)==0)[0]
    pr_unparasitized = len(unparasitized)/Nh

    # mean and variance of parasites per host
    pprh_m = np.mean(paraperhost)
    pprh_v = np.var(paraperhost,ddof=1)

    # correlation of interacting trait pairs
    pntcorr = 0
    if len(traitpairs)>0:
        pntcorr = np.corrcoef(np.transpose(traitpairs))[0,1]

    thg = pd.DataFrame([[Nh,Np,zh,zp,vh,vp,pr_unhosted,pr_unparasitized,pprh_m,pprh_v,pntcorr]], 
        columns=["Nh","Np","zh","zp","vh","vp","pr_unhosted","pr_unparasitized","pprh_m","pprh_v","pntcorr"])
    summary_time_series = summary_time_series.append(thg)

summary_time_series.to_csv("~/gsccs-data/time-series.csv")