from random import gauss
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter
import os

# build another data frame filled with ibd plots

minN = 10.0
maxN = 0.0
minz = 0.0
maxz = 0.0
minv = 10.0
maxv = 0.0

txt = "{time:d}"
time_pts = np.arange(600)
for t in time_pts:    

    fname = "~/gsccs-data/ind-data/indData"+txt.format(time = t).zfill(4)+".csv"
    inds = pd.read_csv(fname)


    # cutting burnt edges of brownie
    res = 4 # the intxn radius
    size = 100 # width/height of geographic region
    gooey_middle = (inds.x>3*res) & (inds.x<size-3*res) & (inds.y>3*res) & (inds.y<size-3*res)
    pretained = sum(gooey_middle)/len(inds.x)

    # discr sp by intxn dist and cutting off outer 3 units in ea direction
    bins = res*np.arange(size/res)

    # binning individuals
    inds.insert(5, "xbin", np.zeros(len(inds.x), dtype=int), True)
    inds.insert(6, "ybin", np.zeros(len(inds.x), dtype=int), True)
    for i in range(len(inds.x)):
        inds.xbin[i] = max(np.where(bins<inds.x[i])[0])
        inds.ybin[i] = max(np.where(bins<inds.y[i])[0])

    # abundances
    Nh_binned = np.zeros((int(size/res)-6,int(size/res)-6))
    Np_binned = np.zeros((int(size/res)-6,int(size/res)-6))

    # mean traits
    zh_binned = np.zeros((int(size/res)-6,int(size/res)-6))
    zp_binned = np.zeros((int(size/res)-6,int(size/res)-6))

    # trait variances
    vh_binned = np.zeros((int(size/res)-6,int(size/res)-6))
    vp_binned = np.zeros((int(size/res)-6,int(size/res)-6))

    xs = np.zeros((int(size/res)-6,int(size/res)-6))
    ys = np.zeros((int(size/res)-6,int(size/res)-6))

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
            if Nh_binned[x,y] > 1:
                vh_binned[x,y] = np.var(inds.z[hs],ddof=1)
            
            if Np_binned[x,y] > 0:
                zp_binned[x,y] = np.mean(inds.z[ps])        
            if Np_binned[x,y] > 1:
                vp_binned[x,y] = np.var(inds.z[ps],ddof=1)

    ref_res = 8*len(bins)-48
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

            xc = int(np.floor(x/8))
            yc = int(np.floor(y/8))
            xs_blurred[x,y] = x
            ys_blurred[x,y] = y

            Nh_blurred[x,y] = Nh_binned[xc,yc]
            zh_blurred[x,y] = zh_binned[xc,yc]
            vh_blurred[x,y] = vh_binned[xc,yc]
            
            Np_blurred[x,y] = Np_binned[xc,yc]
            zp_blurred[x,y] = zp_binned[xc,yc]
            vp_blurred[x,y] = vp_binned[xc,yc]
    
    Nh_blurred = gaussian_filter(Nh_blurred,sigma=4)
    zh_blurred = gaussian_filter(zh_blurred,sigma=4)
    vh_blurred = gaussian_filter(vh_blurred,sigma=4)
    Np_blurred = gaussian_filter(Np_blurred,sigma=4)
    zp_blurred = gaussian_filter(zp_blurred,sigma=4)
    vp_blurred = gaussian_filter(vp_blurred,sigma=4)

    rastdata = {
        'N':np.concatenate((Nh_blurred.flatten(), Np_blurred.flatten())),
        'z':np.concatenate((zh_blurred.flatten(), zp_blurred.flatten())),
        'v':np.concatenate((vh_blurred.flatten(), vp_blurred.flatten())),
        'x':np.concatenate((3.5*res+res*xs_blurred.flatten()/8,3.5*res+res*xs_blurred.flatten()/8)),
        'y':np.concatenate((3.5*res+res*ys_blurred.flatten()/8,3.5*res+res*ys_blurred.flatten()/8))}
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