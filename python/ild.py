#
# computes matrix of interspecific spatial covariance of allele frequencies
#

# note: some rows or columns of cov-ild may be zero as there may be some polymorphic loci
#   with derived type at such low freq's it doesn't appear in the "sample"
#   but this just means the "true" cov for that pair of loci is close to zero, so zero is a good approximation
#   hence, i will keep those zeros so they are reflected in the ild distribution

import os
import numpy as np
import pandas as pd
# os.chdir("../python")
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

# compute ild across whole genome
ild_cor = cf.ild(sstm,"cor")
ild_cov = cf.ild(sstm,"cov")
np.savetxt(fname('ild-cor.csv'), ild_cor, delimiter=",")
np.savetxt(fname('ild-cov.csv'), ild_cov, delimiter=",")
np.savetxt(fname('ild-cor-flat.csv'), ild_cor.flatten(), delimiter=",")
np.savetxt(fname('ild-cov-flat.csv'), ild_cov.flatten(), delimiter=",")

# compute ild only at neutrl loci
hncsl = len(sstm.h.csl_snp)
pncsl = len(sstm.p.csl_snp)
ntl_ild_cor = ild_cor[hncsl:sstm.h.S,pncsl:sstm.p.S]
ntl_ild_cov = ild_cov[hncsl:sstm.h.S,pncsl:sstm.p.S]
np.savetxt(fname('ntl-ild-cor.csv'), ntl_ild_cor, delimiter=",")
np.savetxt(fname('ntl-ild-cov.csv'), ntl_ild_cov, delimiter=",")
np.savetxt(fname('ntl-ild-cor-flat.csv'), ntl_ild_cor.flatten(), delimiter=",")
np.savetxt(fname('ntl-ild-cov-flat.csv'), ntl_ild_cov.flatten(), delimiter=",")

# compute ild only at causal loci
csl_ild_cor = ild_cor[0:hncsl,0:pncsl]
csl_ild_cov = ild_cov[0:hncsl,0:pncsl]
np.savetxt(fname('csl-ild-cor.csv'), csl_ild_cor, delimiter=",")
np.savetxt(fname('csl-ild-cov.csv'), csl_ild_cov, delimiter=",")
np.savetxt(fname('csl-ild-cor-flat.csv'), csl_ild_cor.flatten(), delimiter=",")
np.savetxt(fname('csl-ild-cov-flat.csv'), csl_ild_cov.flatten(), delimiter=",")