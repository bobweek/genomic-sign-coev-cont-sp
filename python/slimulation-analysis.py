import pandas as pd
import numpy as np

# selection parameters
g  = 0.1
sₕ = 0.0
sₚ = 0.01

time_series = pd.DataFrame(columns=[
            "pretained","Nh","Np","pr_unhosted","pr_unparasitized","pprh_m","pprh_v",
            "pntcorr","dsccorr","unocc_h","unocc_p","Nh_m","Np_m","Nh_stdv","Np_stdv",
            "Nh_cv","Np_cv","Ncorr","rh_m","rp_m","rh_stdv","rp_stdv","rh_rp",
            "zh_m", "zh_stdv", "zp_m", "zp_stdv", "zcorr",
            "Dzh_m","Dzp_m","Dzh_stdv","Dzp_stdv","Dzh_Dzp",
            "DBzh_m","DBzp_m","DBzh_stdv","DBzp_stdv","DBzh_DBzp","Nh_DBzp"])

txt = "{time:d}"
time_pts = np.arange(1000,step=10)+1
for t in time_pts:    

    fname = "python/ind-data/indData"+txt.format(time = t).zfill(4)+".csv"
    inds = pd.read_csv(fname)

    hosts = inds.index[inds.spp==1]
    zₕ = inds.z[inds.spp==1]
    zₚ = inds.z[inds.spp==2]

    Bₕ = np.zeros(len(inds.x))
    Bₚ = np.zeros(len(inds.x)) + 1 # any unhosted parasite has unity biotic component of fitness
    traitpairs = []
    paraperhost = []
    hostedparas = np.zeros(len(inds.x))
    for h in hosts:
        paraon = (inds.spp==2) & (inds.x==inds.x[h]) & (inds.y==inds.y[h])    
        paraperhost.append(sum(paraon))
        BH = 1
        for p in inds.index[paraon]:
            hostedparas[p] = 1
            traitpairs.append([zₕ[h],zₚ[p]])
            a = np.exp(-g*(zₕ[h]-zₚ[p])**2/2)
            Bₚ[p] = np.exp(sₚ)*a + (1-a)
            BH *= np.exp(sₕ)*a + (1-a)
        Bₕ[h] = BH

    # cutting burnt edges of brownie
    iota = 5 # the intxn radius
    size = 100 # width/height of geographic region
    gooey_middle = (inds.x>3*iota) & (inds.x<size-3*iota) & (inds.y>3*iota) & (inds.y<size-3*iota)
    pretained = sum(gooey_middle)/len(inds.x)
    hosts = inds.index[(inds.spp==1) & gooey_middle]

    # global census sizes
    Nh = len(hosts)
    Np = len(inds.index[(inds.spp==2) & gooey_middle])

    # proportion of unhosted parasites
    pr_unhosted = (Np - sum(hostedparas[gooey_middle]))/Np

    # proportion of unparasitized hosts
    unparasitized = np.where(np.asarray(paraperhost)==0)[0]
    pr_unparasitized = len(unparasitized)/Nh

    # mean and variance of parasites per host
    pprh_m = np.mean(paraperhost)
    pprh_v = np.var(paraperhost,ddof=1)

    # correlation of interacting trait pairs
    pntcorr = np.corrcoef(np.transpose(traitpairs))[0,1]

    # discr sp by intxn dist and cutting off outer 3 units in ea direction
    bins = iota*np.arange(size/iota)

    # binning individuals
    inds.insert(5, "xbin", np.zeros(len(inds.x), dtype=int), True)
    inds.insert(6, "ybin", np.zeros(len(inds.x), dtype=int), True)
    for i in range(len(inds.x)):
        inds.xbin[i] = max(np.where(bins<inds.x[i])[0])
        inds.ybin[i] = max(np.where(bins<inds.y[i])[0])

    # abundances
    Nh_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))
    Np_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))

    # mean traits
    zh_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))
    zp_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))

    # trait variances
    vh_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))
    vp_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))

    # growth rate
    rh_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))
    rp_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))

    # expected phenotypic response to selection
    Dzh_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))-10
    Dzp_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))-10

    # expected phenotypic response to biotic selection
    DBzh_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))-10
    DBzp_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))-10

    # selection gradients
    bh_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))-10
    bp_binned = np.zeros((int(size/iota)-6,int(size/iota)-6))-10

    for x in range(len(bins)-6):
        for y in range(len(bins)-6):

            hs = (inds.spp==1) & (inds.xbin==x+3) & (inds.ybin==y+3)
            ps = (inds.spp==2) & (inds.xbin==x+3) & (inds.ybin==y+3)

            Nh_binned[x,y] = sum(hs)
            Np_binned[x,y] = sum(ps)

            if Nh_binned[x,y] > 0:
                zh_binned[x,y] = np.mean(inds.z[hs])            
                rh_binned[x,y] = np.mean(np.log(inds.W[hs]))            
            if Nh_binned[x,y] > 1:
                hmat = np.cov(np.transpose(inds[hs][["zh","W"]]))
                Dzh_binned[x,y] = hmat[1,0]
                hBmat = np.cov([np.asarray(zₕ[hs]),Bₕ[hs]])            
                DBzh_binned[x,y] = hBmat[1,0]
                vh_binned[x,y] = hmat[0,0]            
                bh_binned[x,y] = hmat[1,0]/hmat[0,0]
            
            if Np_binned[x,y] > 0:
                zp_binned[x,y] = np.mean(inds.z[ps])        
                rp_binned[x,y] = np.mean(np.log(inds.W[ps]))
            if Np_binned[x,y] > 1:
                pmat = np.cov(np.transpose(inds[ps][["zh","W"]]))
                Dzp_binned[x,y] = pmat[1,0]
                pBmat = np.cov([np.asarray(zₚ[ps]),Bₚ[ps]])
                DBzp_binned[x,y] = pBmat[1,0]            
                vp_binned[x,y] = pmat[0,0]
                bp_binned[x,y] = pmat[1,0]/pmat[0,0]

    # cells occupied by host/parasite resp.
    hocc = Nh_binned.flatten() > 0
    pocc = Np_binned.flatten() > 0

    # cells occupied by several host/parasite resp.
    hsev = Nh_binned.flatten() > 1
    psev = Np_binned.flatten() > 1

    # fraction unoccupied by each spp
    unocc_h = 1-sum(hocc)/((len(bins)-6)**2)
    unocc_p = 1-sum(pocc)/((len(bins)-6)**2)

    # spatial correlation of mean traits across cells
    # note this will result in significantly higher correlations
    # due to increased spatial auto-correlation created by spatial discretization
    dsccorr = np.corrcoef(zh_binned.flatten()[hocc & pocc],zp_binned.flatten()[hocc & pocc])[0,1]

    # mean and var of abund per square for each spp
    Nh_m = np.mean(Nh_binned.flatten())
    Np_m = np.mean(Np_binned.flatten())
    Nh_stdv = np.sqrt(np.var(Nh_binned.flatten(),ddof=1))
    Np_stdv = np.sqrt(np.var(Np_binned.flatten(),ddof=1))
    Nh_cv = Nh_stdv/Nh_m
    Np_cv = Np_stdv/Np_m
    Ncorr = np.corrcoef([Nh_binned.flatten(), Np_binned.flatten()])[0,1]

    # mean and var of local mean trait per square for each spp
    zh_m = np.mean(zh_binned.flatten()[hocc])
    zp_m = np.mean(zp_binned.flatten()[pocc])
    zh_stdv = np.sqrt(np.var(zh_binned.flatten()[hocc],ddof=1))
    zp_stdv = np.sqrt(np.var(zp_binned.flatten()[pocc],ddof=1))
    zcorr = np.corrcoef([Nh_binned.flatten()[hocc & pocc], Np_binned.flatten()[hocc & pocc]])[0,1]

    # mean, stdv and corr of growth rates
    rh_hocc = rh_binned.flatten()[hocc]
    rp_pocc = rp_binned.flatten()[pocc]
    rh_hpocc = rh_binned.flatten()[hocc & pocc]
    rp_hpocc = rp_binned.flatten()[hocc & pocc]
    rh_m = np.mean(rh_hocc)
    rp_m = np.mean(rp_pocc)
    rh_stdv = np.sqrt(np.var(rh_hocc,ddof=1))
    rp_stdv = np.sqrt(np.var(rp_pocc,ddof=1))
    rh_rp = np.corrcoef([rh_hpocc, rp_hpocc])[0,1]

    # mean, var and corr of expected phenotypic change in resp to sel
    Dzh_hocc = Dzh_binned.flatten()[hsev]
    Dzp_pocc = Dzp_binned.flatten()[psev]
    Dzh_hpsev = Dzh_binned.flatten()[hsev & psev]
    Dzp_hpsev = Dzp_binned.flatten()[hsev & psev]
    Dzh_m = np.mean(Dzh_hocc)
    Dzp_m = np.mean(Dzp_pocc)
    Dzh_stdv = np.sqrt(np.var(Dzh_hocc,ddof=1))
    Dzp_stdv = np.sqrt(np.var(Dzp_pocc,ddof=1))
    Dzh_Dzp = np.corrcoef([Dzh_hpsev,Dzp_hpsev])[0,1]

    # mean, var and corr of expected phenotypic change in resp to bio sel
    DBzh_hocc = DBzh_binned.flatten()[hsev]
    DBzp_pocc = DBzp_binned.flatten()[psev]
    DBzh_hpsev = DBzh_binned.flatten()[hsev & psev]
    DBzp_hpsev = DBzp_binned.flatten()[hsev & psev]
    DBzh_m = np.mean(DBzh_hocc)
    DBzp_m = np.mean(DBzp_pocc)
    DBzh_stdv = np.sqrt(np.var(DBzh_hocc,ddof=1))
    DBzp_stdv = np.sqrt(np.var(DBzp_pocc,ddof=1))
    DBzh_DBzp = np.corrcoef([DBzh_hpsev,DBzp_hpsev])[0,1]

    # correlation of host abundance and direction of parasite trait evolution
    Nh_DBzp = np.corrcoef(Nh_binned.flatten()[pocc],DBzp_binned.flatten()[pocc])[0,1]

    # since we set host biotic sel to zero, this should explain parasite trait evolution tracking host density
    hden_ptrt = np.corrcoef(Nh_binned.flatten(),zp_binned.flatten())[0,1]

    # likely won't be helpful unti host experiences coev sel
    bh_hocc = bh_binned.flatten()[hsev]
    bh_m = np.mean(bh_hocc)
    bh_v = np.var(bh_hocc,ddof=1)

    # this could be helpful tho
    bp_pocc = bp_binned.flatten()[psev]
    Np_pocc = Np_binned.flatten()[psev]
    Nh_pocc = Nh_binned.flatten()[psev]
    bp_m = np.mean(bp_pocc)
    bp_v = np.var(bp_pocc,ddof=1)
    bp_Np = np.corrcoef([bp_pocc, Np_pocc])[0,1]
    bp_Nh = np.corrcoef([bp_pocc, Nh_pocc])[0,1]

    # likely won't be helpful until host experiences coev selection
    bh_hpocc = bh_binned.flatten()[hsev & psev]
    bp_hpocc = bp_binned.flatten()[hsev & psev]
    bcorr = np.corrcoef([bh_hpocc, bp_hpocc])[0,1]

    thg = pd.DataFrame([[
        pretained,Nh,Np,pr_unhosted,pr_unparasitized,pprh_m,pprh_v,
        pntcorr,dsccorr,unocc_h,unocc_p,Nh_m,Np_m,Nh_stdv,Np_stdv,
        Nh_cv,Np_cv,Ncorr,rh_m,rp_m,rh_stdv,rp_stdv,rh_rp,
        zh_m, zh_stdv, zp_m, zp_stdv, zcorr,
        Dzh_m,Dzp_m,Dzh_stdv,Dzp_stdv,Dzh_Dzp,
        DBzh_m,DBzp_m,DBzh_stdv,DBzp_stdv,DBzh_DBzp,Nh_DBzp]], 
        columns=[
            "pretained","Nh","Np","pr_unhosted","pr_unparasitized","pprh_m","pprh_v",
            "pntcorr","dsccorr","unocc_h","unocc_p","Nh_m","Np_m","Nh_stdv","Np_stdv",
            "Nh_cv","Np_cv","Ncorr","rh_m","rp_m","rh_stdv","rp_stdv","rh_rp",
            "zh_m", "zh_stdv", "zp_m", "zp_stdv", "zcorr",
            "Dzh_m","Dzp_m","Dzh_stdv","Dzp_stdv","Dzh_Dzp",
            "DBzh_m","DBzp_m","DBzh_stdv","DBzp_stdv","DBzh_DBzp","Nh_DBzp"])
    time_series = time_series.append(thg)

time_series.to_csv("r/time-series.csv")