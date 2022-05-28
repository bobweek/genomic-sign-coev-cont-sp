import pandas as pd
import os
import numpy as np
from sklearn.neighbors import KDTree

# ild = np.load(os.path.expanduser("~/gsccs-data/ild"))

ild = pd.read_csv(os.path.expanduser("~/gsccs-data/iscaf-cov.csv"), header=None)
ild = ild.to_numpy()

ild.shape

csl_ild = pd.read_csv(os.path.expanduser("~/gsccs-data/csl-iscaf-cov.csv"), header=None)
csl_ild = csl_ild.to_numpy()

# remove host neutral loci within 5kbp of host causal loci

h_ntl_snps = pd.read_csv(os.path.expanduser("~/gsccs-data/h-ntl-snps.csv"), header=None).to_numpy().flatten()
h_csl_snps = pd.read_csv(os.path.expanduser("~/gsccs-data/h-csl-snps.csv"), header=None).to_numpy().flatten()

hntlKD = KDTree(h_ntl_snps.reshape(-1,1))
hntlnrbs = np.zeros(0)
for i in np.arange(len(h_csl_snps)):
    nbri = hntlKD.query_radius(h_csl_snps[i].reshape(1,-1),r=5e3)[0]
    hntlnrbs = np.concatenate((hntlnrbs,nbri))

# remove para neutral loci within 5kbp of para causal loci

p_ntl_snps = pd.read_csv(os.path.expanduser("~/gsccs-data/p-ntl-snps.csv"), header=None).to_numpy().flatten()
p_csl_snps = pd.read_csv(os.path.expanduser("~/gsccs-data/p-csl-snps.csv"), header=None).to_numpy().flatten()

pntlKD = KDTree(p_ntl_snps.reshape(-1,1))
pntlnrbs = np.zeros(0)
for i in np.arange(len(p_csl_snps)):
    nbri = pntlKD.query_radius(p_csl_snps[i].reshape(1,-1),r=5e3)[0]
    pntlnrbs = np.concatenate((pntlnrbs,nbri))

ild = ild.flatten()
csl_ild = csl_ild.flatten()

ild = ild[ild!=0]
csl_ild = csl_ild[csl_ild!=0]

thrsh = 1e-4

upper = sum(csl_ild > np.quantile(ild,1-thrsh))
lower = sum(csl_ild < np.quantile(ild,thrsh))
ttl_frac = (upper+lower)/len(csl_ild)