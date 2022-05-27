import numpy as np
import allel
import pandas as pd
from scipy import spatial
from dataclasses import dataclass

# model parameters dataclass
@dataclass
class Pars:
    iota: float     # spatial interspp intxn dist
    # eventually ill add in other parameters
    def __iter__(self):
        return iter((self.iota))

# position and genotype information dataclass for a given species
@dataclass
class Species:
    csl_gt: allel.GenotypeArray # genotype array of the species
    ntl_gt: allel.GenotypeArray # genotype array of the species
    loc: np.ndarray         # geographic locations of individuals
    csl_snp: np.ndarray     # genomic positions of causal snps
    ntl_snp: np.ndarray     # genomic positions of neutral snps
    snp: np.ndarray         # genomic positions of all
    S: int                  # number polymorphic loci
    csl_frq: np.ndarray     # locus-by-ind matrix of causal freqs
    ntl_frq: np.ndarray     # locus-by-ind matrix of neutral freqs
    frq: np.ndarray         # locus-by-ind matrix of all freqs (csl staked on top)
    trt: np.ndarray         # trait values of individuals
    def __iter__(self):
        return iter((self.gt, self.loc, self.csl_snp, self.ntl_snp, self.trt))

# position and genotype information dataclass for whole system
@dataclass
class System:
    h: Species  # host species dataclass
    p: Species  # para species dataclass
    pars: Pars  # model parameters
    def __iter__(self):
        return iter((self.h, self.p, self.pars))

def loadUp(hmet,pmet,hga_csl,pga_csl,hga_ntl,pga_ntl,
    hcsl_snp,pcsl_snp,hntl_snp,pntl_snp,params,hcsl_frq,pcsl_frq,hntl_frq,pntl_frq):

    # read in geo locations and trait values for each sample
    h_metadat = pd.read_csv(hmet, delimiter = ",")
    p_metadat = pd.read_csv(pmet, delimiter = ",")

    # pull out geo locations and trait values
    h_loc = np.array(h_metadat[['x','y']])
    p_loc = np.array(p_metadat[['x','y']])
    h_trt = np.array(h_metadat['z'])
    p_trt = np.array(p_metadat['z'])

    # read in snp positions
    h_csl_snp = np.loadtxt(hcsl_snp, delimiter=",")
    p_csl_snp = np.loadtxt(pcsl_snp, delimiter=",")
    
    # read in snp positions
    h_ntl_snp = np.loadtxt(hntl_snp, delimiter=",")
    p_ntl_snp = np.loadtxt(pntl_snp, delimiter=",")

    # combine snp positions (with causal in front)
    h_snp = np.concatenate((h_csl_snp,h_ntl_snp))
    p_snp = np.concatenate((p_csl_snp,p_ntl_snp))

    # pull out genotype arrays of causal mutations
    h_csl_gt = allel.GenotypeArray(np.load(hga_csl))
    p_csl_gt = allel.GenotypeArray(np.load(pga_csl))

    # pull out genotype arrays of neutral mutations
    h_ntl_gt = allel.GenotypeArray(np.load(hga_ntl))
    p_ntl_gt = allel.GenotypeArray(np.load(pga_ntl))

    # compute number polymorphic loci
    hS = h_csl_gt.shape[0] + h_ntl_gt.shape[0]
    pS = p_csl_gt.shape[0] + p_ntl_gt.shape[0]

    # get locus-by-individual matrix of allele freqs
    hcslfrq = np.load(hcsl_frq)
    pcslfrq = np.load(pcsl_frq)
    hntlfrq = np.load(hntl_frq)
    pntlfrq = np.load(pntl_frq)
    hfrq = np.concatenate((hcslfrq, hntlfrq))
    pfrq = np.concatenate((pcslfrq, pntlfrq))

    # grab model parameters
    prs = pd.read_csv(params, delimiter = ",")

    # combine genotype arrays, geo locations,
    # traits and snp positions
    # into one dataclass to rule them all
    h = Species(csl_gt=h_csl_gt, ntl_gt=h_ntl_gt, loc=h_loc, trt=h_trt, 
        csl_snp=h_csl_snp, ntl_snp=h_ntl_snp, snp=h_snp, S=hS, csl_frq=hcslfrq, ntl_frq=hntlfrq, frq=hfrq)
    p = Species(csl_gt=p_csl_gt, ntl_gt=p_ntl_gt, loc=p_loc, trt=p_trt, 
        csl_snp=p_csl_snp, ntl_snp=p_ntl_snp, snp=p_snp, S=pS, csl_frq=pcslfrq, ntl_frq=pntlfrq, frq=pfrq)
    pars = Pars(iota=prs["ι"][0])
    sys = System(h=h,p=p,pars=pars)

    return sys

# interspecific spatial covariance of allele frequencies
def ild(sstm,co):

    # first build kd-tree for host locations
    # for every parasite, select a random host within radius iota
    # this forms an "interacting pair"
    # compute the spatial covariance of allele freqs across pairs
    # repeating this process will redraw interacting pairs, thus allowing precision assessment
    # will be v slow tho

    h, p, pars = sstm

    h_tree = spatial.KDTree(h.loc)
    
    h_h = []
    for ploc in p.loc:
        h_nbrs = h_tree.query_ball_point(ploc,3*pars.iota)
        if len(h_nbrs)>0:
            h_h.append(np.random.choice(h_nbrs))
        else:
            # indicate when no host neighbors
            h_h.append(-1)
    
    # remove -1's from choisen host indeces and
    # remove parasites that didn't find a host
    h_h = np.array(h_h)
    intx_indxs = np.where(h_h>-1)[0]
    hs = h_h[intx_indxs]
    ps = np.arange(len(p.trt))[intx_indxs]

    if co=="cov":
        
        # spatial cov allele freq
        ild = np.zeros((h.S,p.S))    
        for i in np.arange(h.S):
            for j in np.arange(p.S):
                h_freqs = h.frq[i,hs]
                p_freqs = p.frq[j,ps]
                ild[i,j] = np.cov(h_freqs,p_freqs)[0,1]

    else:
        
        # spatial cor allele freq
        ild = np.zeros((h.S,p.S))    
        for i in np.arange(h.S):
            for j in np.arange(p.S):
                h_freqs = h.frq[i,hs]
                p_freqs = p.frq[j,ps]
                ild[i,j] = np.corrcoef(h_freqs,p_freqs)[0,1]
                if np.isnan(ild[i,j]):
                    ild[i,j] = 0

    return ild
