import os
import numpy as np
import classes_fcts as cf

def iname(dpth,nme):
    pth = dpth+nme    
    return os.path.expanduser(pth)

def fname(nme):
    pth = "~/gsccs-data/"+nme    
    return os.path.expanduser(pth)

def ild(dpth,r):

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
    params = 'slim-pars.csv'
    hcsl_frq = fname('hfrq-causal.npy')
    pcsl_frq = fname('pfrq-causal.npy')
    hntl_frq = fname('hfrq-neutrl.npy')
    pntl_frq = fname('pfrq-neutrl.npy')

    sstm = cf.loadUp(hmet,pmet,hga_csl,pga_csl,hga_ntl,pga_ntl,
        hcsl_snp,pcsl_snp,hntl_snp,pntl_snp,params,hcsl_frq,pcsl_frq,hntl_frq,pntl_frq)

    # compute iscaf across whole genome (only for cov-based iscaf)
    iscaf = cf.iscaf(sstm,"cov")
    np.savetxt(iname(dpth,'iscaf%i.csv'%r), iscaf, delimiter=",")
    np.savetxt(iname(dpth,'iscaf-flat%i.csv'%r), iscaf.flatten(), delimiter=",")

    # compute iscaf only at causal loci
    hncsl = len(sstm.h.csl_snp)
    pncsl = len(sstm.p.csl_snp)
    csl_iscaf = iscaf[0:hncsl,0:pncsl]
    np.savetxt(iname(dpth,'csl-iscaf%i.csv'%r), csl_iscaf, delimiter=",")
    np.savetxt(iname(dpth,'csl-iscaf-flat%i.csv'%r), csl_iscaf.flatten(), delimiter=",")