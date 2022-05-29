import burnin as bp
import genotypearrays as ga
import ild
import numpy as np
import pandas as pd
import subprocess
import multiprocessing as mp
import os

#
# this script first populates the subfolder structure that organizes
# parameter combinations and replicates per parameter combo
# with parameter files that input to msprime and slim
# then it will run the makeILD function in parallel across all
# parameter combos and replicates, but just ncpu at a time
# in turn, the makeILD function populates the subfolder structure
# with initial tree seqs generated by msprime, with final tree seqs
# generated by slim, separate genotype arrays for causal and neutral
# mutations, and, finally, causal, neutral and combined ild matrices
#

# btw, this creates the subfolder structure used w/in folders sxs and Lxs
# for i in {0..8}; do mkdir $i; cd $i; for j in {0..4}; do mkdir $j; done; cd ../; done

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
mus = cprs["mu"][1:3]
ks = cprs["k"]

def makeILD(j,r,wch):
    
    # path to msprime and slim model parameters in superfolder
    msprpth = os.path.expanduser("~/gsccs-data/replicates/"+wch+"/%i"%j+"/")
    
    # replicate specific subfolder
    datapth = msprpth+"%i"%r+"/"

    # do the burn-in
    bp.burnin(msprpth,datapth)

    # run the slimulation calling parameter file in specified subfolder
    slmcmd = """slim -d "parfname='"""+msprpth+"""slim-pars.csv'"  -d "trefldr='"""+datapth+"""'" hp.slim"""
    subprocess.run(slmcmd, shell=True)

    # add neutral mutations and export genotype arrays
    ga.genotypeArrays(datapth)

    # compute ild matrices and save to specified subfolder  
    ild.ild(datapth,msprpth)

# the thing that manages parallel running of the things
pool = mp.Pool(5)

# # iterate across combinations of host-para biotic selection
# j=0
# for sh in ss:
#     for sp in ss:
        
#         superfolder = os.path.expanduser("~/gsccs-data/replicates/sxs/%i"%j)
        
#         # swap out parameters
#         sprs["sₕ"] = sh
#         sprs["sₚ"] = sp
#         bprs.to_csv(superfolder+"/bparams.csv", index=False)
#         sprs.to_csv(superfolder+"/slim-pars.csv", index=False)

#         [pool.apply_async(makeILD, args=(j,r,"sxs")) for r in np.arange(reps)]
        
#         j+=1

# # iterate across combinations of para biotic selection and mutation rate
# j=0
# for mu in mus:
#     for sp in ss:
        
#         superfolder = os.path.expanduser("~/gsccs-data/replicates/Lxs/%i"%j)

#         # swap out parameters
#         wchk = np.where(cprs["mu"] == mu)[0][0]
#         k = ks[wchk]
#         bprs["mu"] = mu
#         sprs["sₚ"] = sp
#         bprs["k"] = k
#         sprs["κₚ"] = k
#         bprs.to_csv(superfolder+"/bparams.csv", index=False)
#         sprs.to_csv(superfolder+"/slim-pars.csv", index=False)
        
#         [pool.apply_async(makeILD, args=(j,r,"Lxs")) for r in np.arange(reps)]

#         j+=1

# run sub loop to get remaining Lxs'
# iterate across combinations of para biotic selection and mutation rate
j=3
for mu in mus:
    for sp in ss:
        
        superfolder = os.path.expanduser("~/gsccs-data/replicates/Lxs/%i"%j)

        # swap out parameters
        k = np.float64(cprs["k"][cprs["mu"] == mu])
        bprs["mu"] = mu
        sprs["sₚ"] = sp
        bprs["k"] = k
        sprs["κₚ"] = k
        bprs.to_csv(superfolder+"/bparams.csv", index=False)
        sprs.to_csv(superfolder+"/slim-pars.csv", index=False)
        
        [pool.apply_async(makeILD, args=(j,r,"Lxs")) for r in np.arange(reps)]

        j+=1

pool.close()
pool.join()