import numpy as np
import classes_fcts as cf

def ild(dpth,ppth):

    params = ppth+'slim-pars.csv'

    hmet = dpth+'host.txt'
    pmet = dpth+'para.txt'
    hga_csl = dpth+'hga-causal.npy'
    pga_csl = dpth+'pga-causal.npy'
    hga_ntl = dpth+'hga-neutrl.npy'
    pga_ntl = dpth+'pga-neutrl.npy'
    hcsl_snp = dpth+'h-csl-snps.csv'
    pcsl_snp = dpth+'p-csl-snps.csv'
    hntl_snp = dpth+'h-ntl-snps.csv'
    pntl_snp = dpth+'p-ntl-snps.csv'    
    hcsl_frq = dpth+'hfrq-causal.npy'
    pcsl_frq = dpth+'pfrq-causal.npy'
    hntl_frq = dpth+'hfrq-neutrl.npy'
    pntl_frq = dpth+'pfrq-neutrl.npy'

    sstm = cf.loadUp(hmet,pmet,hga_csl,pga_csl,hga_ntl,pga_ntl,
        hcsl_snp,pcsl_snp,hntl_snp,pntl_snp,params,hcsl_frq,pcsl_frq,hntl_frq,pntl_frq)

    # compute ild across whole genome (only for cov-based ild)
    ild = cf.ild(sstm,"cov")
    np.savetxt(dpth+'ild.csv', ild, delimiter=",")
    np.savetxt(dpth+'ild-flat.csv', ild.flatten(), delimiter=",")

    # compute ild only at causal loci
    hncsl = len(sstm.h.csl_snp)
    pncsl = len(sstm.p.csl_snp)
    csl_ild = ild[0:hncsl,0:pncsl]
    np.savetxt(dpth+'csl-ild.csv', csl_ild, delimiter=",")
    np.savetxt(dpth+'csl-ild-flat.csv', csl_ild.flatten(), delimiter=",")

    # compute ild only at neutral loci
    ntl_ild = ild[hncsl:sstm.h.S,pncsl:sstm.p.S]
    np.savetxt(dpth+'ntl-ild.csv', ntl_ild, delimiter=",")
    np.savetxt(dpth+'ntl-ild-flat.csv', ntl_ild.flatten(), delimiter=",")