import numpy as np
import allel
import pandas as pd
from scipy import stats
from scipy import spatial
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from sklearn.neighbors import KernelDensity
from dataclasses import dataclass
import matplotlib.pyplot as plt 

# model parameters dataclass
@dataclass
class Pars:
    SI: float       # spatial interspp intxn dist
    hc: float       # multiplicative cost on host fitness for each infection
    pc: float       # multiplicative cost on parasite fitness for failing to infect
    pb: float       # multiplicative benefit on parasite fitness for successful infection
    minpr: float    # minimum probability of infection
    maxpr: float    # maximum probability of infection
    lmbda: float    # width of spatial niche
    ch: float       # strength of host competition
    cp: float       # strength of para competition
    def __iter__(self):
        return iter((self.SI, self.hc, self.pc, self.pb, self.minpr, 
            self.maxpr, self.lmbda, self.ch, self.cp))

# position and genotype information dataclass for a given species
@dataclass
class Species:
    gt: allel.GenotypeArray # genotype array of the species
    pos: np.ndarray         # geo positions of individuals
    trt: np.ndarray         # trait values of individuals
    causL: int              # genomic position of causal locus
    def __iter__(self):
        return iter((self.gt, self.pos, self.trt, self.causL))

# position and genotype information dataclass for whole system
@dataclass
class System:
    h: Species  # host species dataclass
    p: Species  # para species dataclass
    def __iter__(self):
        return iter((self.h, self.p))

# same as Species, but with discretized spatial positions
@dataclass
class discSpecies:
    p: np.ndarray       # allele frequencies in each cell
    G: np.ndarray       # genetic variance in each cell
    N: np.ndarray       # abundances in each cell
    trt: list           # trait values of individuals in each cell
    causL: int          # genomic position of causal locus
    S: int              # number of variant loci being tracked
    def __iter__(self):
        return iter((self.f_c, self.f_n, self.N, self.trt, self.causL, self.S))

# same as System, but with discretized spatial positions
@dataclass
class discSystem:
    h: discSpecies  # host species dataclass
    p: discSpecies  # para species dataclass
    def __iter__(self):
        return iter((self.h, self.p))

def loadUp(hmet,pmet,causL,hvcf,pvcf):
    # read in positions and trait values for each sample
    h_metadat = pd.read_csv(hmet, delimiter = "\t")
    p_metadat = pd.read_csv(pmet, delimiter = "\t")

    # read in genomic positions of causal loci
    causL_dat = pd.read_csv(causL, delimiter = "\t")
    h_causL = causL_dat["h"][0]
    p_causL = causL_dat["p"][0]

    # read in vcf files
    h_callset = allel.read_vcf(hvcf)
    p_callset = allel.read_vcf(pvcf)

    # pull out genotype arrays
    h_gt = allel.GenotypeArray(h_callset['calldata/GT'])
    p_gt = allel.GenotypeArray(p_callset['calldata/GT'])

    # pull out positions and trait values
    h_pos = np.array(h_metadat[['x','y']])
    p_pos = np.array(p_metadat[['x','y']])
    h_trt = np.array(h_metadat['z'])
    p_trt = np.array(p_metadat['z'])

    # combine genotype arrays, geo positions,
    # traits and positions of causal loci
    # into one dataclass to rule them all
    h = Species(gt=h_gt, pos=h_pos, trt=h_trt, causL=h_causL)
    p = Species(gt=p_gt, pos=p_pos, trt=p_trt, causL=p_causL)
    sys = System(h=h,p=p)

    return sys

# estimates selection coefficients using random samples from each species
def glbl_selCoeffs(sys,pars,h_pr,p_pr):

    # unpack data
    # note: gt & causL are not used
    h, p = sys
    h_gt, h_pos, h_trt, h_causL = h
    p_gt, p_pos, p_trt, p_causL = p

    # unpack pars
    SI, hc, pc, pb, minpr, maxpr, lmbda, ch, cp = pars

    # do the calc
    h_N = len(h_trt)
    p_N = len(p_trt)
    h_n = int(np.ceil(h_pr*h_N))
    p_n = int(np.ceil(p_pr*p_N))
    h_smpl = np.random.choice(np.arange(h_N),h_n,replace=False)
    p_smpl = np.random.choice(np.arange(p_N),p_n,replace=False)
    h_tree = spatial.KDTree(h_pos)
    p_tree = spatial.KDTree(p_pos)
    h_fits = []
    for h in h_smpl:
        p_nbrs = p_tree.query_ball_point(h_pos[h],SI)
        h_fit = []
        # accumulate interspp intxns
        for p in p_nbrs:        
            h_inf = h_tree.query_ball_point(p_pos[p],SI)
            p_prb = 1/len(h_inf)        # pr that p attempts to infect focal host ind
            if h_trt[h] == p_trt[p]:    # case of matching genotype
                h_fit.append(hc*maxpr*p_prb+(1-maxpr*p_prb)) # avg of hc and no effect
            else:                       # case of mis-matched genotype
                h_fit.append(hc*minpr*p_prb+(1-minpr*p_prb))
        if len(p_nbrs) == 0:
            h_fit = 1.0                 # if no parasites around, no effect
        h_fit = np.mean(h_fit)
        # accumulate intraspp intxns
        h_nbrs = h_tree.query_ball_point(h_pos[h],3*2*lmbda)
        h_comp = len(h_nbrs)
        # h_comp = 0
        # for hnb in h_nbrs:
        #     d = spatial.distance.euclidean(h_pos[h],h_pos[hnb])
        #     h_comp += U*U*stats.norm.pdf(d,0,np.sqrt(2*lmbda))    # spatial niche overlap..... (eqn SM.80)
        h_fit *= np.exp(-ch*h_comp)                                 # spatial competition effect (eqn SM.82)
        h_fits.append(h_fit)
    h_fits = np.array(h_fits)
    h_0 = np.where(h_trt[h_smpl]==0)[0] # indices of genotypes
    h_1 = np.where(h_trt[h_smpl]==1)[0]
    h_2 = np.where(h_trt[h_smpl]==2)[0]
    h_W00 = np.mean(h_fits[h_0]) # mean fit of ea genotype
    h_W01 = np.mean(h_fits[h_1])
    h_W11 = np.mean(h_fits[h_2])
    h_N0 = len(h_0) # allele counts
    h_N1 = len(h_1)
    h_N2 = len(h_2)
    h_W0 = (2*h_W00*h_N0 + h_W01*h_N1)/(2*h_N0+h_N1) # mean fit of ea allele
    h_W1 = (2*h_W11*h_N2 + h_W01*h_N1)/(2*h_N2+h_N1)
    h_s = h_W1 - h_W0 # host selection coefficient

    p_fits = []
    for p in p_smpl:
        # accumulate interspp intxns
        h_nbrs = h_tree.query_ball_point(p_pos[p],SI)
        p_fit = []
        for h in h_nbrs:
            if p_trt[p] == h_trt[h]:
                p_fit.append(pb*maxpr+pc*(1-maxpr)) # avg of benefit from infecting and cost of not infecting
            else:
                p_fit.append(pb*minpr+pc*(1-minpr))
        if len(h_nbrs) == 0:
            p_fit = pc                              # cost of not having anyone to infect
        p_fit = np.mean(p_fit)
        # accumulate intraspp intxns
        p_nbrs = p_tree.query_ball_point(p_pos[p],3*2*lmbda)
        p_comp = len(p_nbrs)
        # p_comp = 0
        # for pnb in p_nbrs:
        #     d = spatial.distance.euclidean(p_pos[p],p_pos[pnb])
        #     p_comp += U*U*stats.norm.pdf(d,0,np.sqrt(2*lmbda)) # spatial niche overlap..... (eqn SM.80)
        p_fit *= np.exp(-cp*p_comp)                 # spatial competition effect (eqn SM.82)
        p_fits.append(p_fit)
    p_fits = np.array(p_fits)
    p_0 = np.where(p_trt[p_smpl]==0)[0] # indices of genotypes
    p_1 = np.where(p_trt[p_smpl]==1)[0]
    p_2 = np.where(p_trt[p_smpl]==2)[0]
    p_W00 = np.mean(p_fits[p_0]) # mean fit of ea genotype
    p_W01 = np.mean(p_fits[p_1])
    p_W11 = np.mean(p_fits[p_2])
    p_N0 = len(p_0) # allele counts
    p_N1 = len(p_1)
    p_N2 = len(p_2)
    p_W0 = (2*p_W00*p_N0 + p_W01*p_N1)/(2*p_N0+p_N1) # mean fit of ea allele
    p_W1 = (2*p_W11*p_N2 + p_W01*p_N1)/(2*p_N2+p_N1)
    p_s = p_W1 - p_W0 # parasite selection coefficient

    return h_s, p_s

def dropSmolFreqs(sys, p_min):

    # unpack data
    h, p = sys
    h_gt, h_pos, h_trt, h_causL = h
    p_gt, p_pos, p_trt, p_causL = p

    h_glb_ac = h_gt.count_alleles()
    h_ply1 = h_glb_ac[:,0]/(h_glb_ac[:,0]+h_glb_ac[:,1]) > p_min # ancestl allele freq at least f_min
    h_ply2 = h_glb_ac[:,1]/(h_glb_ac[:,0]+h_glb_ac[:,1]) > p_min # derived allele freq at least f_min
    h_poly = [a and b for a, b in zip(h_ply1.tolist(), h_ply2.tolist())] # both anc and der at least f_min
    h_poly[h_causL] = False # make sure throw away causal locus (since the actual information is saved in trait metadata)
    h_gt = h_gt[h_poly,:,:]

    p_glb_ac = p_gt.count_alleles()
    p_ply1 = p_glb_ac[:,0]/(p_glb_ac[:,0]+p_glb_ac[:,1]) > p_min # ancestl allele freq at least f_min
    p_ply2 = p_glb_ac[:,1]/(p_glb_ac[:,0]+p_glb_ac[:,1]) > p_min # derived allele freq at least f_min
    p_poly = [a and b for a, b in zip(p_ply1.tolist(), p_ply2.tolist())]
    p_poly[p_causL] = False # make sure throw away causal locus
    p_gt = p_gt[p_poly,:,:]

    # recalculate genomic positions of causal loci
    h_causL = sum(h_poly[0:h_causL])
    p_causL = sum(p_poly[0:p_causL])

    # pack it up
    h = Species(gt=h_gt, pos=h_pos, trt=h_trt, causL=h_causL)
    p = Species(gt=p_gt, pos=p_pos, trt=p_trt, causL=p_causL)
    sys = System(h=h,p=p)

    return sys

def discSpace(sys, width, height, res):

    # unpack data
    h, p = sys
    h_gt, h_pos, h_trt, h_causL = h
    p_gt, p_pos, p_trt, p_causL = p

    # bin space
    hor_bins = width*np.arange(res)/res
    ver_bins = height*np.arange(res)/res

    # in each cell, compute abundance N, frequencies of causal alleles, and gather allele counts at neutral loci
    h_S = h_gt.n_variants
    p_S = p_gt.n_variants
    hbn_N = np.zeros((res,res)) # matrix containing abundance of inds in ea cell
    pbn_N = np.zeros((res,res)) # same, but for parasite
    hbn_p = np.zeros((h_S+1,res,res)) # contains frequencies of neutral alleles in corr bin
    pbn_p = np.zeros((p_S+1,res,res)) # the additional locus will be for splicing in causal freq

    hbn_caus = np.zeros((res,res)) # matrix containing freq of causal allele in ea cell
    pbn_caus = np.zeros((res,res)) # same, but for parasite
    hbn_trt = [] # a 'matrix' whose entries are lists of traits in corresponding cell
    pbn_trt = [] # same, but for parasite
    hbn_gt = [] # a 'matrix' where entries are list of host allele cnt matrices in corresponding cell
    pbn_gt = [] # same, but for parasite
    for X in hor_bins:
        hbn_trt.append([])
        pbn_trt.append([])
        hbn_gt.append([])
        pbn_gt.append([])
        h_tlist = []
        h_glist = []
        p_tlist = []
        p_glist = []
        for Y in ver_bins:
            indX = np.where(hor_bins==X)[0][0]
            indY = np.where(ver_bins==Y)[0][0]
            h_in_bnX1 = X < h_pos[:,0]
            h_in_bnX2 = h_pos[:,0] < X+width/res
            h_in_bnX = [a and b for a, b in zip(h_in_bnX1.tolist(), h_in_bnX2.tolist())]
            h_in_bnY1 = Y < h_pos[:,1]
            h_in_bnY2 = h_pos[:,1] < Y+height/res
            h_in_bnY = [a and b for a, b in zip(h_in_bnY1.tolist(), h_in_bnY2.tolist())]
            h_mmbrs = [a and b for a, b in zip(h_in_bnX, h_in_bnY)]
            h_tlist = h_trt[h_mmbrs]
            h_glist = h_gt[:,h_mmbrs]
            hbn_N[indX,indY] = len(h_tlist)
            hbn_trt[indX].append(h_tlist)
            hbn_gt[indX].append(h_glist)
            if len(h_tlist)>0:
                hbn_caus[indX,indY] = sum(h_tlist)/(2*len(h_tlist))
            else:
                hbn_caus[indX,indY] = np.nan
            p_in_bnX1 = X < p_pos[:,0]
            p_in_bnX2 = p_pos[:,0] < X+width/res
            p_in_bnX = [a and b for a, b in zip(p_in_bnX1.tolist(), p_in_bnX2.tolist())]
            p_in_bnY1 = Y < p_pos[:,1]
            p_in_bnY2 = p_pos[:,1] < Y+height/res
            p_in_bnY = [a and b for a, b in zip(p_in_bnY1.tolist(), p_in_bnY2.tolist())]
            p_mmbrs = [a and b for a, b in zip(p_in_bnX, p_in_bnY)]
            p_tlist = p_trt[p_mmbrs]
            p_glist = p_gt[:,p_mmbrs]
            pbn_N[indX,indY] = len(p_tlist)
            pbn_trt[indX].append(p_tlist)
            pbn_gt[indX].append(p_glist)
            if len(p_tlist)>0:      # check that this cell is occupied, else place nan
                pbn_caus[indX,indY] = sum(p_tlist)/(2*len(p_tlist))
            else:
                pbn_caus[indX,indY] = np.nan

    # compute allele frequencies of neutral sites for each cell
    for i in np.arange(res):
        for j in np.arange(res):
            if hbn_N[i][j] > 0: # makes sure there's inds in that cell (else places a nan)

                h_ac = hbn_gt[i][j].count_alleles()
                h_f = h_ac[:,1]/(h_ac[:,0]+h_ac[:,1]) # computes the freqs
                
                # splice in causal allele freq
                hbn_p[0:h_causL,i,j] = h_f[0:h_causL]
                hbn_p[h_causL,i,j] = hbn_caus[i,j]
                hbn_p[h_causL+1:h_S+1,i,j] = h_f[h_causL:h_S]
            else:
                hbn_p[:,i,j] = np.nan
            
            if pbn_N[i][j] > 0:
                p_ac = pbn_gt[i][j].count_alleles()
                p_f = p_ac[:,1]/(p_ac[:,0]+p_ac[:,1])

                # splice in causal allele freq
                pbn_p[0:p_causL,i,j] = p_f[0:p_causL]
                pbn_p[p_causL,i,j] = pbn_caus[i,j]
                pbn_p[p_causL+1:p_S+1,i,j] = p_f[p_causL:p_S]
            else:
                pbn_p[:,i,j] = np.nan
    
    # compute genetic variance at each locus in each bin
    h_G = hbn_p*(1-hbn_p)
    p_G = pbn_p*(1-pbn_p)

    # pack it up
    h = discSpecies(p=hbn_p, G=h_G, N=hbn_N, trt=hbn_trt, causL=h_causL, S=(h_S+1))
    p = discSpecies(p=pbn_p, G=p_G, N=pbn_N, trt=pbn_trt, causL=p_causL, S=(p_S+1))
    dis_sys = discSystem(h=h,p=p)

    return dis_sys

def Fst(dis_sys):
    h, p = dis_sys

    h_pbar = np.zeros(h.S)
    h_Gbar = np.zeros(h.S)
    for i in np.arange(h.S):
        h_pbar[i] = np.average(h.p[:,:][i], weights=h.N)
        Gbars = h.p[:,:][i]*(1-h.p[:,:][i])
        h_Gbar[i] = np.average(Gbars, weights=h.N)

    p_pbar = np.zeros(p.S)
    p_Gbar = np.zeros(p.S)
    for i in np.arange(p.S):
        p_pbar[i] = np.average(p.p[:,:][i], weights=p.N)
        Gbars = p.p[:,:][i]*(1-p.p[:,:][i])
        p_Gbar[i] = np.average(Gbars, weights=p.N)

    h_G = h_pbar*(1-h_pbar)
    h_Fst = 1-h_Gbar/h_G

    p_G = p_pbar*(1-p_pbar)
    p_Fst = 1-p_Gbar/p_G

    return h_Fst, p_Fst

def Fst_plts(h_Fst,p_Fst,dis_sys):

    h, p = dis_sys

    # "manhattan" plot of Fst for host
    plt.figure()
    plt.scatter(x=np.arange(h.S),y=h_Fst)
    plt.scatter(x=h.causL,y=h_Fst[h.causL],c="red")
    plt.savefig('h_Fst.png')

    # "manhattan" plot of Fst for parasite
    plt.figure()
    plt.scatter(x=np.arange(p.S),y=p_Fst)
    plt.scatter(x=p.causL,y=p_Fst[p.causL],c="red")
    plt.savefig('p_Fst.png')

# return indices of values in the k-th percentile 
def singleSppFilter(stat, k, disSp):

    filtered = np.where(stat > np.percentile(stat,k))[0]
    if len(np.where(filtered==disSp.causL)[0]) > 0:
        causL = np.where(filtered==disSp.causL)[0][0]
    else:
        causL = np.nan
        print("warning: causal locus lost")
    p = disSp.p[filtered]
    G = disSp.G[filtered]
    S = len(filtered)

    # pack it up
    newSp = discSpecies(p=p,G=G,N=disSp.N,trt=disSp.trt,causL=causL,S=S)

    return newSp

def discaf(dis_sys):

    h, p = dis_sys

    # spatially corr allele freq
    iscaf = np.zeros((h.S,p.S))
    
    # find where both species are present
    h_presence = h.N.flatten() > 0
    p_presence = p.N.flatten() > 0
    presence = [a and b for a, b in zip(h_presence.tolist(), p_presence.tolist())]
    
    # compute correlations across cells where both are present
    for i in np.arange(h.S):
        for j in np.arange(p.S):
            h_freqs = h.p[i].flatten()[presence]
            p_freqs = p.p[j].flatten()[presence]
            iscaf[i,j] = np.corrcoef(h_freqs,p_freqs)[0,1]

    return iscaf

def causal_search(iscaf, dis_sys):

    h, p = dis_sys

    maxcorr = np.where(iscaf == max(iscaf.flatten()))
    if h.causL==maxcorr[0][0] and p.causL==maxcorr[1][0]:
        print("causal loci found!")
    else:
        print("causal loci NOT found!")


def iscaf_plts(iscaf,csl_iscaf):

    # plot matrix of iscafs
    plt.figure()
    matplt = plt.matshow(iscaf)
    plt.colorbar(matplt)
    plt.savefig('iscaf_matrix.png')
    plt.close()

    # make density plot of iscaf
    # iscaf_kde = stats.gaussian_kde(iscaf.flatten())
    # xrang = np.arange(0,1,0.01)
    # plt.figure()
    # plt.plot(xrang,iscaf_kde.pdf(xrang))
    # if not np.isnan(csl_iscaf):
    #     plt.scatter(x=csl_iscaf,y=iscaf_kde.pdf(csl_iscaf),c="red")
    # plt.savefig('iscaf_kde.png')
    # plt.close()

    # make histogram of iscaf    
    plt.figure()    
    plt.hist(iscaf.flatten(),bins=100,alpha=0.5)
    plt.scatter(x=csl_iscaf,y=0,c="red")
    plt.savefig('iscaf_hist.png')
    plt.close()

def csl_mnhtn_plts(iscaf,dis_sys):

    h, p = dis_sys

    # "manhattan" plot of host iscaf with para causal locus
    plt.figure()
    plt.scatter(x=np.arange(h.S),y=iscaf[:,p.causL])
    plt.scatter(x=h.causL,y=iscaf[h.causL,p.causL],c="red")
    plt.savefig('h_iscaf_at_causal_p.png')
    plt.close()

    # "manhattan" plot of para iscaf with host causal locus
    plt.figure()
    plt.scatter(x=np.arange(p.S),y=iscaf[h.causL,:])
    plt.scatter(x=p.causL,y=iscaf[h.causL,p.causL],c="red")
    plt.savefig('p_iscaf_at_causal_h.png')
    plt.close()
    
def interp_p(s,l,res):

    grid_x, grid_y = np.mgrid[0:res:1, 0:res:1]
    notnans = np.where(np.isnan(s.p[l])==False)
    pts = np.column_stack((notnans[0],notnans[1]))
    values = s.p[l][notnans]
    intrp = griddata(pts, values, (grid_x,grid_y), method='cubic')

    return intrp

def plot_p(dis_sys,res):
    h, p = dis_sys

    # interpolate missing values for causal loci
    h_srf = [interp_p(h,h.causL,res)]
    p_srf = [interp_p(p,p.causL,res)]

    # three random loci for each species
    cands = np.delete(np.arange(h.S),h.causL)
    loci = np.random.choice(cands,3,replace=False)
    for l in loci:
        # interpolate missing values
        h_srf.append(interp_p(h,l,res))

    cands = np.delete(np.arange(p.S),p.causL)
    loci = np.random.choice(cands,3,replace=False)
    for l in loci:
        p_srf.append(interp_p(p,l,res))

    # plot host allele freq surfs
    plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.suptitle('Host Allele Freq Surfaces')
    l = 0
    for ax in axes.flat:
        im = ax.matshow(h_srf[l], vmin=0, vmax=1)
        l += 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig('h_freq_surf.png')
    plt.close()

    # plot para allele freq surfs
    plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.suptitle('Parasite Allele Freq Surfaces')
    l = 0
    for ax in axes.flat:
        im = ax.matshow(p_srf[l], vmin=0, vmax=1)
        l += 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig('p_freq_surf.png')
    plt.close()

    #
    # now do the same thing, but smooth each surface so its easier to look at
    #

    # plot host allele freq surfs
    plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.suptitle('Host Allele Freq Surfaces')
    l = 0
    for ax in axes.flat:
        blurred = gaussian_filter(h_srf[l], sigma=0.5)
        im = ax.matshow(blurred, vmin=0, vmax=1)
        l += 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig('h_freq_surf_blr.png')
    plt.close()

    # para host allele freq surfs
    plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.suptitle('Parasite Allele Freq Surfaces')
    l = 0
    for ax in axes.flat:
        blurred = gaussian_filter(p_srf[l], sigma=0.5)
        im = ax.matshow(blurred, vmin=0, vmax=1)
        l += 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig('p_freq_surf_blr.png')
    plt.close()

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)