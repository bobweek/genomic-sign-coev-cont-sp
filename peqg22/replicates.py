import burnin as bp
import genotypearrays as ga
import ild
import numpy as np
import pandas as pd
import subprocess
import multiprocessing as mp

# # in case you need to reload module
import importlib
importlib.reload(ild) # example use

# number of replicates per parameter combo
reps = 5

# burn-in parameters
bprs = pd.read_csv("bparams.csv", sep=",")

# slim parameters
sprs = pd.read_csv("slim-pars.csv", sep=",")

# parameter combos to iterate over
cprs = pd.read_csv("parcombos.csv", sep=",")

# pull out combos
ss = cprs["s"]
mus = cprs["mu"]

def makeILD(j,r,wch):

    # subfolder for specified parameter combo (replicates w same pars get pooled together)
    datapth = "~/gsccs-data/replicates/"+wch+"/%i"%j+"/"

    # do the burn-in
    bp.burnin(datapth+"bparams.csv")

    # run the slimulation calling parameter file in specified subfolder
    slmcmd = """slim -d "parfname='"""+datapth+"""slim-pars.csv'" hp.slim"""
    subprocess.run(slmcmd, shell=True)

    # add neutral mutations and export genotype arrays
    ga.genotypeArrays()

    # compute ild matrices and save to specified subfolder    
    ild.ild(datapth,r)

# the thing that manages parallel running of the things
pool = mp.Pool()

# iterate across combinations of host-para biotic selection
j=0
for sh in ss:
    for sp in ss:

        j+=1
        
        # swap out parameters
        sprs["sₕ"] = sh
        sprs["sₚ"] = sp
        bprs.to_csv("~/gsccs-data/replicates/sxs/%i"%j+"/bparams.csv", index=False)
        sprs.to_csv("~/gsccs-data/replicates/sxs/%i"%j+"/slim-pars.csv", index=False)

        [pool.apply_async(makeILD, args=(j,r,"sxs")) for r in np.arange(reps)]
        
        # save par combo to specified subfolder
        sprs.to_csv("~/gsccs-data/replicates/sxs/%i"%j+"/pars.csv", index=False)        

# iterate across combinations of para biotic selection and mutation rate
j=0
for mu in mus:
    for sp in ss:

        j+=1

        # swap out parameters
        bprs["mu"] = mu
        sprs["sₚ"] = sp
        bprs.to_csv("~/gsccs-data/replicates/Lxs/%i"%j+"/bparams.csv", index=False)
        sprs.to_csv("~/gsccs-data/replicates/Lxs/%i"%j+"/slim-pars.csv", index=False)
        
        [pool.apply_async(makeILD, args=(j,r,"Lxs")) for r in np.arange(reps)]

        # save par combo to specified subfolder
        sprs.to_csv("~/gsccs-data/replicates/Lxs/%i"%j+"/pars.csv", index=False)        

# restore default mutation rate
bprs["mu"] = 5.0e-13
bprs.to_csv("bparams.csv", index=False)