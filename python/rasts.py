from random import gauss
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter
import os

# selection parameters
g  = 0.1
sₕ = 0.0
sₚ = 0.01

# build another data frame filled with ibd plots

minN = 10.0
maxN = 0.0
minz = 0.0
maxz = 0.0
minv = 10.0
maxv = 0.0

txt = "{time:d}"
time_pts = np.arange(1001)
for t in time_pts:    

    fname = "~/gsccs-data/ind-data/indData"+txt.format(time = t).zfill(4)+".csv"
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

    xs = np.zeros((int(size/iota)-6,int(size/iota)-6))
    ys = np.zeros((int(size/iota)-6,int(size/iota)-6))

    for x in range(len(bins)-6):
        for y in range(len(bins)-6):

            xs[x,y] = x
            ys[x,y] = y

            hs = (inds.spp==1) & (inds.xbin==x+3) & (inds.ybin==y+3)
            ps = (inds.spp==2) & (inds.xbin==x+3) & (inds.ybin==y+3)

            Nh_binned[x,y] = sum(hs)
            Np_binned[x,y] = sum(ps)

            if Nh_binned[x,y] > 0:
                zh_binned[x,y] = np.mean(inds.z[hs])            
                rh_binned[x,y] = np.mean(np.log(inds.W[hs]))            
            if Nh_binned[x,y] > 1:
                hmat = np.cov(np.transpose(inds[hs][["z","W"]]))
                Dzh_binned[x,y] = hmat[1,0]
                # hBmat = np.cov([np.asarray(zₕ[hs]),Bₕ[hs]])            
                # DBzh_binned[x,y] = hBmat[1,0]
                vh_binned[x,y] = hmat[0,0]            
            
            if Np_binned[x,y] > 0:
                zp_binned[x,y] = np.mean(inds.z[ps])        
                rp_binned[x,y] = np.mean(np.log(inds.W[ps]))
            if Np_binned[x,y] > 1:
                pmat = np.cov(np.transpose(inds[ps][["z","W"]]))
                Dzp_binned[x,y] = pmat[1,0]
                # pBmat = np.cov([np.asarray(zₚ[ps]),Bₚ[ps]])
                # DBzp_binned[x,y] = pBmat[1,0]            
                vp_binned[x,y] = pmat[0,0]

    ref_res = 4*len(bins)-24
    xs_blurred = np.zeros((ref_res,ref_res))
    ys_blurred = np.zeros((ref_res,ref_res))
    Nh_blurred = np.zeros((ref_res,ref_res))
    zh_blurred = np.zeros((ref_res,ref_res))
    vh_blurred = np.zeros((ref_res,ref_res))
    Np_blurred = np.zeros((ref_res,ref_res))
    zp_blurred = np.zeros((ref_res,ref_res))
    vp_blurred = np.zeros((ref_res,ref_res))
    for x in range(ref_res):
        for y in range(ref_res):

            xc = int(np.floor(x/4))
            yc = int(np.floor(y/4))
            xs_blurred[x,y] = x
            ys_blurred[x,y] = y

            Nh_blurred[x,y] = Nh_binned[xc,yc]
            zh_blurred[x,y] = zh_binned[xc,yc]
            vh_blurred[x,y] = vh_binned[xc,yc]
            
            Np_blurred[x,y] = Np_binned[xc,yc]
            zp_blurred[x,y] = zp_binned[xc,yc]
            vp_blurred[x,y] = vp_binned[xc,yc]
    
    Nh_blurred = gaussian_filter(Nh_blurred,sigma=1)
    zh_blurred = gaussian_filter(zh_blurred,sigma=1)
    vh_blurred = gaussian_filter(vh_blurred,sigma=1)
    Np_blurred = gaussian_filter(Np_blurred,sigma=1)
    zp_blurred = gaussian_filter(zp_blurred,sigma=1)
    vp_blurred = gaussian_filter(vp_blurred,sigma=1)

    rastdata = {
        'N':np.concatenate((Nh_blurred.flatten(), Np_blurred.flatten())),
        'z':np.concatenate((zh_blurred.flatten(), zp_blurred.flatten())),
        'v':np.concatenate((vh_blurred.flatten(), vp_blurred.flatten())),
        'x':np.concatenate((3.5*iota+iota*xs_blurred.flatten()/4,3.5*iota+iota*xs_blurred.flatten()/4)),
        'y':np.concatenate((3.5*iota+iota*ys_blurred.flatten()/4,3.5*iota+iota*ys_blurred.flatten()/4))}
    rastdf = pd.DataFrame(rastdata)
    
    minN = min(minN,min(rastdf.N))
    maxN = max(maxN,max(rastdf.N))
    minz = min(minz,min(rastdf.z))
    maxz = max(maxz,max(rastdf.z))
    minv = min(minv,min(rastdf.v))
    maxv = max(maxv,max(rastdf.v))

    rname = "~/gsccs-data/rast-data/rast"+txt.format(time = t).zfill(4)+".csv"
    rastdf.to_csv(rname)

minmaxdata = {
    'minN':np.asarray([minN]),
    'maxN':np.asarray([maxN]),
    'minz':np.asarray([minz]),
    'maxz':np.asarray([maxz]),
    'minv':np.asarray([minv]),
    'maxv':np.asarray([maxv])
}
minmaxdf = pd.DataFrame(minmaxdata)
minmaxdf.to_csv("~/gsccs-data/rast-data/minmax.csv")

print("DONE!")

duration = 1  # seconds
freq = 440  # Hz
os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))