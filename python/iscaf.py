#
# computes matrix of interspecific spatial covariance of allele frequencies
#

# note: some rows or columns of cov-iscaf may be zero as there may be some polymorphic loci
#   with derived type at such low freq's it doesn't appear in the "sample"
#   but this just means the "true" cov for that pair of loci is close to zero, so zero is a good approximation
#   hence, i will keep those zeros so they are reflected in the ild distribution

import os
import numpy as np
import pandas as pd
os.chdir("python")
import classes_fcts as cf

# # in case you need to reload module
# import importlib
# importlib.reload(cf) # example use

def fname(nme):
    pth = '~/gsccs-data/'+nme    
    return os.path.expanduser(pth)

hmet = fname('host.txt')
pmet = fname('para.txt')
hga_csl = fname('hga-causal.npy')
pga_csl = fname('pga-causal.npy')
hga_ntl = fname('hga-neutrl.npy')
pga_ntl = fname('pga-neutrl.npy')
hcsl_snp = fname('h-csl-snps.csv')
pcsl_snp = fname('p-csl-snps.csv')
hntl_snp = fname('h-ntl-snps.csv')
pntl_snp = fname('p-ntl-snps.csv')
params = fname('params.csv')
hcsl_frq = fname('hfrq-causal.npy')
pcsl_frq = fname('pfrq-causal.npy')
hntl_frq = fname('hfrq-neutrl.npy')
pntl_frq = fname('pfrq-neutrl.npy')

prs = pd.read_csv(params, delimiter = ",")

sstm = cf.loadUp(hmet,pmet,hga_csl,pga_csl,hga_ntl,pga_ntl,
    hcsl_snp,pcsl_snp,hntl_snp,pntl_snp,params,hcsl_frq,pcsl_frq,hntl_frq,pntl_frq)

# compute iscaf across whole genome
iscaf_cor = cf.iscaf(sstm,"cor")
iscaf_cov = cf.iscaf(sstm,"cov")
np.savetxt(fname('iscaf-cor.csv'), iscaf_cor, delimiter=",")
np.savetxt(fname('iscaf-cov.csv'), iscaf_cov, delimiter=",")
np.savetxt(fname('iscaf-cor-flat.csv'), iscaf_cor.flatten(), delimiter=",")
np.savetxt(fname('iscaf-cov-flat.csv'), iscaf_cov.flatten(), delimiter=",")

# compute iscaf only at causal loci
hncsl = len(sstm.h.csl_snp)
pncsl = len(sstm.p.csl_snp)
csl_iscaf_cor = iscaf_cor[0:hncsl,0:pncsl]
csl_iscaf_cov = iscaf_cov[0:hncsl,0:pncsl]
np.savetxt(fname('csl-iscaf-cor.csv'), csl_iscaf_cor, delimiter=",")
np.savetxt(fname('csl-iscaf-cov.csv'), csl_iscaf_cov, delimiter=",")
np.savetxt(fname('csl-iscaf-cor-flat.csv'), csl_iscaf_cor.flatten(), delimiter=",")
np.savetxt(fname('csl-iscaf-cov-flat.csv'), csl_iscaf_cov.flatten(), delimiter=",")

# import networkx as nx

# import matplotlib.pyplot as plt

# import sys
# import warnings
# warnings.filterwarnings('ignore')

# print('Python Version : '+sys.version)
# print('NetworkX version : '+nx.__version__)