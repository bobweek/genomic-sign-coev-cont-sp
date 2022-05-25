import os
import pandas as pd
os.chdir("python")
import classes_fcts as cf

# in case you need to reload module
# import importlib
# importlib.reload(cf) # example use


hmet = os.path.expanduser('~/gsccs-data/host.txt')
pmet = os.path.expanduser('~/gsccs-data/para.txt')
hga_csl = os.path.expanduser('~/gsccs-data/hga-causal.npy')
pga_csl = os.path.expanduser('~/gsccs-data/pga-causal.npy')
hga_ntl = os.path.expanduser('~/gsccs-data/hga-neutrl.npy')
pga_ntl = os.path.expanduser('~/gsccs-data/pga-neutrl.npy')
hcsl_snp = os.path.expanduser('~/gsccs-data/h_csl_snps.csv')
pcsl_snp = os.path.expanduser('~/gsccs-data/p_csl_snps.csv')
hntl_snp = os.path.expanduser('~/gsccs-data/h_ntl_snps.csv')
pntl_snp = os.path.expanduser('~/gsccs-data/p_ntl_snps.csv')
params = os.path.expanduser('~/gsccs-data/params.csv')
hcsl_frq = os.path.expanduser('~/gsccs-data/hfrq-causal.npy')
pcsl_frq = os.path.expanduser('~/gsccs-data/pfrq-causal.npy')
hntl_frq = os.path.expanduser('~/gsccs-data/hfrq-neutrl.npy')
pntl_frq = os.path.expanduser('~/gsccs-data/pfrq-neutrl.npy')

prs = pd.read_csv(params, delimiter = ",")

sys = cf.loadUp(hmet,pmet,hga_csl,pga_csl,hga_ntl,pga_ntl,
    hcsl_snp,pcsl_snp,hntl_snp,pntl_snp,params,hcsl_frq,pcsl_frq,hntl_frq,pntl_frq)

iscaf = cf.iscaf(sys)

iscaf.shape