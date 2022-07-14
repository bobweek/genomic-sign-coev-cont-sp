import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

# this needs to be saved somewhere else
reps = 5

threshold = 1e-2

# parameter combos to iterate over
cprs = pd.read_csv("parcombos.csv", sep=",")

# pull out combos
ss = cprs["s"]

def perCsl(j,r,wch,thrsh):

    dpth = os.path.expanduser("~/gsccs-data/replicates/"+wch+"/%i"%j+"/%i"%r+"/")

    # ild = pd.read_csv(dpth+"ild.csv", header=None).to_numpy()
    # csl_ild = pd.read_csv(dpth+"csl-ild.csv", header=None).to_numpy()
    ild = np.load(dpth+"ild.csv.npy")
    csl_ild = np.load(dpth+"csl-ild.npy")

    # remove host neutral loci within 5kbp of host causal loci

    h_ntl_snps = pd.read_csv(dpth+"h-ntl-snps.csv", header=None).to_numpy().flatten()
    h_csl_snps = pd.read_csv(dpth+"h-csl-snps.csv", header=None).to_numpy().flatten()

    hntlKD = KDTree(h_ntl_snps.reshape(-1,1))
    hntlnrbs = np.zeros(0)
    for i in np.arange(len(h_csl_snps)):
        nbri = hntlKD.query_radius(h_csl_snps[i].reshape(1,-1),r=5e3)[0]
        hntlnrbs = np.concatenate((hntlnrbs,nbri))

    ild = np.delete(ild,[int(nbr) for nbr in hntlnrbs],0)

    # remove para neutral loci within 5kbp of para causal loci

    p_ntl_snps = pd.read_csv(dpth+"p-ntl-snps.csv", header=None).to_numpy().flatten()
    p_csl_snps = pd.read_csv(dpth+"p-csl-snps.csv", header=None).to_numpy().flatten()

    pntlKD = KDTree(p_ntl_snps.reshape(-1,1))
    pntlnrbs = np.zeros(0)
    for i in np.arange(len(p_csl_snps)):
        nbri = pntlKD.query_radius(p_csl_snps[i].reshape(1,-1),r=5e3)[0]
        pntlnrbs = np.concatenate((pntlnrbs,nbri))

    ild = np.delete(ild,[int(nbr) for nbr in pntlnrbs],1)

    ild = ild.flatten()
    csl_ild = csl_ild.flatten()

    ild = ild[ild!=0]
    csl_ild = csl_ild[csl_ild!=0]

    thrsh = 1e-2

    upper = sum(csl_ild > np.quantile(ild,1-thrsh))
    lower = sum(csl_ild < np.quantile(ild,thrsh))
    ttl_frac = (upper+lower)/len(csl_ild)

    return ttl_frac

def getL(j,r,wch,S):
    dpth = os.path.expanduser("~/gsccs-data/replicates/"+wch+"/%i"%j+"/%i"%r+"/")
    p_csl_snps = pd.read_csv(dpth+S+"-csl-snps.csv", header=None).to_numpy().flatten()
    return len(p_csl_snps)

# arrays containing results used for heatmaps
sfrac = np.zeros((3,3))
mufrac = np.zeros((3,3))
sss = np.zeros((2,3,3))
Ls = np.zeros((1+reps,3,3))

# populate heatmap arrays
j=0
for k in np.arange(3):
    for l in np.arange(3):
        
        sfrs = [perCsl(j,r,"sxs",threshold) for r in np.arange(reps)]
        sfrac[k,l] = np.mean(sfrs)
        sss[0,k,:] = ss[k] # sh
        sss[1,:,l] = ss[l] # sp
        
        mufrs = [perCsl(j,r,"Lxs",threshold) for r in np.arange(reps)]
        mufrac[k,l] = np.mean(mufrs)
        Ls[0,k,l] = ss[k]
        Ls[1+np.arange(reps),k,l] = [getL(j,r,"Lxs","p") for r in np.arange(reps)]        

        j += 1

superfolder = os.path.expanduser("~/gsccs-data/replicates/")
np.savetxt(superfolder+"sxs/cslfrac.csv",sfrac.flatten())
np.savetxt(superfolder+"Lxs/cslfrac.csv",mufrac.flatten())
for i in np.arange(2):
    np.savetxt(superfolder+"sxs/ss%i.csv"%(i+1),sss[i,:,:].flatten())
np.savetxt(superfolder+"Lxs/s.csv",Ls[0,:,:].flatten())
for r in np.arange(reps):
    np.savetxt(superfolder+"Lxs/L%i.csv"%(r+1),Ls[(r+1),:,:].flatten())

# # handy snippet for populating subfolder structure with dummy data
# import shutil
# wchs = ["sxs","Lxs"]
# for j in np.arange(9):
#     for wch in wchs:
#         superfolder = os.path.expanduser("~/gsccs-data/replicates/"+wch+"/%i"%j+"/")
#         sprfldrs = [superfolder+"%i"%r+"/" for r in np.arange(reps)]

#         dtfldr = os.path.expanduser("~/gsccs-data/")
#         thgs = ["h-csl-snps.csv", "h-ntl-snps.csv", "p-csl-snps.csv", "p-ntl-snps.csv",
#             "iscaf-cov.csv", "csl-iscaf-cov.csv"]
#         for fldr in sprfldrs:
#             for thg in thgs:
#                 shutil.copy(dtfldr+thg,fldr+thg)