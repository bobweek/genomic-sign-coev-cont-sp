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

    # compute ild across whole genome
    ild_cor = cf.ild(sstm,"cor")
    ild_cov = cf.ild(sstm,"cov")
    np.savetxt(dpth+'ild-cor.csv', ild_cor, delimiter=",")
    np.savetxt(dpth+'ild-cov.csv', ild_cov, delimiter=",")
    np.savetxt(dpth+'ild-cor-flat.csv', ild_cor.flatten(), delimiter=",")
    np.savetxt(dpth+'ild-cov-flat.csv', ild_cov.flatten(), delimiter=",")

    # compute ild only at neutrl loci
    hncsl = len(sstm.h.csl_snp)
    pncsl = len(sstm.p.csl_snp)
    ntl_ild_cor = ild_cor[hncsl:sstm.h.S,pncsl:sstm.p.S]
    ntl_ild_cov = ild_cov[hncsl:sstm.h.S,pncsl:sstm.p.S]
    np.savetxt(dpth+'ntl-ild-cor.csv', ntl_ild_cor, delimiter=",")
    np.savetxt(dpth+'ntl-ild-cov.csv', ntl_ild_cov, delimiter=",")
    np.savetxt(dpth+'ntl-ild-cor-flat.csv', ntl_ild_cor.flatten(), delimiter=",")
    np.savetxt(dpth+'ntl-ild-cov-flat.csv', ntl_ild_cov.flatten(), delimiter=",")

    # compute ild only at causal loci
    csl_ild_cor = ild_cor[0:hncsl,0:pncsl]
    csl_ild_cov = ild_cov[0:hncsl,0:pncsl]
    np.savetxt(dpth+'csl-ild-cor.csv', csl_ild_cor, delimiter=",")
    np.savetxt(dpth+'csl-ild-cov.csv', csl_ild_cov, delimiter=",")
    np.savetxt(dpth+'csl-ild-cor-flat.csv', csl_ild_cor.flatten(), delimiter=",")
    np.savetxt(dpth+'csl-ild-cov-flat.csv', csl_ild_cov.flatten(), delimiter=",")