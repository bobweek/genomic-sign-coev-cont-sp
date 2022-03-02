from os import unlink
import time # for timing computations
import allel
import pyslim # for pulling model parameters from metadata
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import mam_classes_fcts as cf
# from sklearn.neighbors import KernelDensity
# from sklearn.model_selection import GridSearchCV
from scipy.stats import t, normaltest, kstest, anderson, shapiro

# for reloading custom module while developing
import importlib
# importlib.reload(cf) # example use

#
# load in data
#

# yields type(sys) = System
sys = cf.loadUp('host.mam.txt','parasite.mam.txt','causL.mam.txt',
                'h_snps.csv','p_snps.csv','host.mam.npy','para.mam.npy')

#
# drop low frequency loci (prevents large spurious correlations when we discretize)
#

# yields type(sys) = System
sys = cf.dropSmolFreqs(sys, p_min=0.2)
np.savetxt("h_midfreq_snps.csv", sys.h.snp, delimiter=",")
np.savetxt("p_midfreq_snps.csv", sys.p.snp, delimiter=",")

#
# identify independent loci in each species to form null distribution(s)
#

hgt_sqsh = sys.h.gt[:,:,0]+sys.h.gt[:,:,1]
hnlnkd = allel.locate_unlinked(hgt_sqsh)

pgt_sqsh = sys.p.gt[:,:,0]+sys.p.gt[:,:,1]
pnlnkd = allel.locate_unlinked(pgt_sqsh)

# splice in flag to ignore causal loci
h_unlinked = np.arange(sys.h.gt.n_variants+1)[np.concatenate((hnlnkd[:sys.h.causL],np.full(1,False),hnlnkd[sys.h.causL:]))]
p_unlinked = np.arange(sys.p.gt.n_variants+1)[np.concatenate((pnlnkd[:sys.p.causL],np.full(1,False),pnlnkd[sys.p.causL:]))]

# NEED TO REMOVE LOCI THAT ARE "TOO CLOSE" TO CAUSAL LOCI

#
# discretize space into res x res grid of cells
#

# spatial resolution
res = 32

# spatial extent
width = 100 # eventually read in from file
height = 100 # eventually read in from file

# discretization, yields type(dis_sys) = discSystem
# dis_sys = cf.discSpace(sys, width, height, res)
# cf.saveDisSys(dis_sys,"discr_sp_data")
# del dis_sys
dis_sys = cf.loadDisSys("discr_sp_data/dis_sys")
h, p = dis_sys

# prevent unnecessary use of memory
del sys

#
# plot allele freq surfaces of causal loci and three randomly selected neutral loci
#

cf.plot_p(dis_sys,res)

#
# filter loci using single spp Fst patterns (hold off for now)
#

# # how does host-para coev affect Fst?
# h_Fst, p_Fst = cf.Fst(dis_sys)

# # make "manhattan" plots
# cf.Fst_plts(h_Fst,p_Fst,dis_sys)

# # coev appears to significantly lower host Fst
# # coev appears to slightly increase para Fst
# # this seems to agree with the surface plots
# # does this pattern generalize?
# # for what region of parameter space?
# # what other patterns occur?

# # filter out host loci in highest 90% of Fst scores
# dis_sys.h = cf.singleSppFilter(-h_Fst, 90, h)

# # filter out para loci in lowest 90% of Fst scores
# dis_sys.p = cf.singleSppFilter(p_Fst, 90, p)

# # update species specific objects
# h, p = dis_sys

#
# compute interspecific spatial correlations of allele frequencies
#

# compute
# time_start = time.clock_gettime(0)
# iscaf = cf.discaf(dis_sys)
# time_elapsed = (time.clock_gettime(0) - time_start)
# np.savetxt("ild.csv", iscaf, delimiter=",")
iscaf = np.loadtxt("ild.csv", delimiter=",")

# pull out iscaf across unlinked loci
iscaf_unlinked = iscaf[h_unlinked,:][:,p_unlinked]

iscaf[h.causL,p.causL] # the causal corr

# THOT: iscaf defines a bipartite network...

# # fit a t-distr
# t_stats = []
# for l_h in np.arange(h.S):
#     for l_p in np.arange(p.S):
#         r = iscaf[l_h,l_p]
#         t_stt = r*np.sqrt((2.5*res-2)/(1-r**2))
#         t_stats.append(t_stt)

# # see if computed t-statistics fit the student t-distr
# # using kolmogorov-smirnov test
# kstest(t_stats,"t",args=(2*res-2,))

# # visual check if t is good fit (NEED TO LOOK AT TAILS!)
# plt.figure()
# xrang = np.arange(min(t_stats),max(t_stats),0.001)
# plt.scatter(x=xrang,y=t.pdf(xrang,2.5*res-2),c="red")
# plt.scatter(x=np.mean(t_stats),y=0,c="black")
# plt.hist(t_stats, bins=100, alpha=0.5, density=True, stacked=True)
# plt.savefig('t_hist_clr.png')
# plt.close()

# are the causal loci correlated more than any other pair?
# cf.causal_search(iscaf, dis_sys)

# make manhattan plots for causal loci
# cf.csl_mnhtn_plts(abs(iscaf),dis_sys)

#
# distribution of iscaf across independent (w/in spp) loci
#

# make histogram of iscaf and superimpose causal iscaf
# cf.iscaf_plts(iscaf_unlinked,iscaf[h.causL,p.causL])

#
# does student's t-distr apply?
#

# # check if allele frequencies are normally distributed across locations
# # using three different tests
# h_shaps = []
# h_andrs = []
# h_nrmlt = []
# for l in h_unlinked:
#     h_shaps.append(shapiro(h.p[l,:,:].flatten()).pvalue)
#     h_andrs.append(anderson(h.p[l,:,:].flatten()).significance_level)
#     h_nrmlt.append(normaltest(h.p[l,:,:].flatten()).pvalue)
# p_shaps = []
# p_andrs = []
# p_nrmlt = []
# for l in p_unlinked:
#     p_shaps.append(shapiro(p.p[l,:,:].flatten()).pvalue)
#     p_andrs.append(anderson(p.p[l,:,:].flatten()).significance_level)
#     p_nrmlt.append(normaltest(p.p[l,:,:].flatten()).pvalue)

# # just using shapiro-wilks test for now
# h_normies = np.where(np.array(h_shaps) < 0.05)[0]
# p_normies = np.where(np.array(p_shaps) < 0.05)[0]
h_normies = h_unlinked
p_normies = p_unlinked

# proportion of loci ignoring
# worth finding a transform to try and normalize (logit)?
# 1-len(h_normies)/len(h_shaps)
# 1-len(p_normies)/len(p_shaps)

#
# removing intraspecific spatial correlations of allele freqs
# so we can apply students t
#

# estimate covariance function
# this assumes allele freq surfs are sampled from same distribution
# so we can average spatial power spectra across loci to obtain cov

from scipy import signal, linalg

# should all unlinked neutral loci have same expected allele freq?

# this creates res-by-res matrix of distances as function of
# horiz/vert displacement along disc-grid (max = sqrt(2)*res/2)
# since boundaries are periodic, each row/col is unimodal
# this matrix will be important for computing the correlation-adjusted
# cross-correlation matrix (eqn 9 of jendoubi & strimmer 2019)
dists = cf.makeDistMat(res)
d = np.unique(dists)

# analog of dists but for freq domain
freqs = cf.makeFreqMat(res)
f = np.unique(freqs)

#
# here we 
#   1) estimate per locus cov fcts
#   2) fit exponential curves to them to get
#       2.1) per locus spatial scale
#       2.2) per locus marginal var
#   3) average cov fcts across neutral loci to obtain model for "removing space" from candidate loci
#

# return vector of cov/psd's in order of incr dist/freq
hcov_fcts = np.zeros((len(h_normies),len(d)))
pcov_fcts = np.zeros((len(p_normies),len(d)))
hpsd_fcts = np.zeros((len(h_normies),len(f)))
ppsd_fcts = np.zeros((len(p_normies),len(f)))

# host power spectra and cov fcts
hfft = np.zeros((len(h_normies),res,res//2+1), dtype=complex)
hps = np.zeros((len(h_normies),res,res//2+1)) # ps = power spectrum
hcovs = np.zeros((len(h_normies),res,res)) # cov fct per locus
Vhs = np.zeros(len(h_normies))  # marginal variance per locus
xihs = np.zeros(len(h_normies)) # spatial scale per locus
k = 0
for l in h_normies:    
    hps[k,:,:], hcovs[k,:,:] = cf.fft_cov(h.p[l,:,:],h.p[l,:,:])
    for dst in np.arange(len(d)):         # est cov fct
        coords = np.where(dists==d[dst])
        n = len(coords[0])
        h_sum = 0
        for m in np.arange(n):
            i = coords[0][m]
            j = coords[1][m]
            h_sum += hcovs[k,i,j]
        hcov_fcts[k,dst] = h_sum/n
    for frq in np.arange(len(f)):
        coords = np.where(freqs==f[frq])
        n = len(coords[0])
        h_sum = 0
        for m in np.arange(n):
            i = coords[0][m]
            j = coords[1][m]
            h_sum += hps[k,i,j]
        hpsd_fcts[k,frq] = h_sum/n
    # tail end of cov fct tends be garbage
    # so only look up to first negative point
    # to estimate the spatial scale
    neg = min(np.where(hcov_fcts[k,:]<0)[0])
    Vhs[k] = hcovs[k,0,0]
    y = np.log(hcov_fcts[k,1:neg]/Vhs[k])
    xihs[k] = -np.mean(d[1:neg]/y)    # fit exponential cov fct
    k += 1
hpsd = np.mean(hps,axis=0)
hcov = np.fft.irfft2(hpsd) # this is the averaged cov fct for "de-spacing"

# parasite power spectra and cov fcts
pps = np.zeros((len(p_normies),res,res//2+1))
pcovs = np.zeros((len(p_normies),res,res))
Vps = np.zeros(len(p_normies))  # marginal variance per locus
xips = np.zeros(len(p_normies)) # spatial scale per locus
k = 0
for l in p_normies:
    pps[k,:,:], pcovs[k,:,:] = cf.fft_cov(p.p[l,:,:],p.p[l,:,:])
    for dst in np.arange(len(d)):         # est cov fct
        coords = np.where(dists==d[dst])
        n = len(coords[0])
        p_sum = 0
        for m in np.arange(n):
            i = coords[0][m]
            j = coords[1][m]
            p_sum += pcovs[k,i,j]
        pcov_fcts[k,dst] = p_sum/n
    for frq in np.arange(len(f)): # replace with a build_fct function
        coords = np.where(freqs==f[frq])
        n = len(coords[0])
        p_sum = 0
        for m in np.arange(n):
            i = coords[0][m]
            j = coords[1][m]
            p_sum += pps[k,i,j]
        ppsd_fcts[k,frq] = p_sum/n
    neg = min(np.where(pcov_fcts[k,:]<0)[0])
    Vps[k] = pcovs[k,0,0]
    y = np.log(pcov_fcts[k,1:neg]/Vps[k])
    xips[k] = -np.mean(d[1:neg]/y)    # fit exponential cov fct
    k += 1
ppsd = np.mean(pps,axis=0)
pcov = np.fft.irfft2(ppsd) # this is the averaged cov fct for "de-spacing"

# correlation-adjusted cross-correlation matrix K
#   - first estimate cov fct C from unlinked neutral loci
#   - then build cov matrix S for spatial locations
#   - from S compute correlation matrix P = inv(V^0.5).S.inv(V^0.5)
#   - where V^0.5 is diag matrix with values sqrt(C(0)) on diag
#   - for each interspp pair of loci, approx cross-cov fct Chp
#   - build cross-cov matrix Shp and obtain Php = inv(Vh^0.5).Shp.inv(Vp^0.5)
#   - with corr matrices Ph, Pp compute K = inv(Ph^0.5).Php.inv(Pp^0.5)
# this will be a res**2-by-res**2 matrix!!


hcsl_ps, hcsl_cov = cf.fft_cov(h.p[h.causL,:,:],h.p[h.causL,:,:])
pcsl_ps, pcsl_cov = cf.fft_cov(p.p[p.causL,:,:],p.p[p.causL,:,:])
hcov_fct = np.zeros(len(d))
pcov_fct = np.zeros(len(d))
hpsd_fct = np.zeros(len(f))
ppsd_fct = np.zeros(len(f))
hcsl_cov_fct = np.zeros(len(d))
pcsl_cov_fct = np.zeros(len(d))
hcsl_psd_fct = np.zeros(len(f))
pcsl_psd_fct = np.zeros(len(f))
for dst in np.arange(len(d)):
    coords = np.where(dists==d[dst])
    n = len(coords[0])
    h_sum = 0
    p_sum = 0
    hcsl_sum = 0
    pcsl_sum = 0
    for m in np.arange(n):
        i = coords[0][m]
        j = coords[1][m]
        h_sum += hcov[i,j]
        p_sum += pcov[i,j]
        hcsl_sum += hcsl_cov[i,j]
        pcsl_sum += pcsl_cov[i,j]
    hcov_fct[dst] = h_sum/n
    pcov_fct[dst] = p_sum/n
    hcsl_cov_fct[dst] = hcsl_sum/n
    pcsl_cov_fct[dst] = pcsl_sum/n
for frq in np.arange(len(f)): # replace with a build_fct function
    coords = np.where(freqs==f[frq])
    n = len(coords[0])
    h_sum = 0
    p_sum = 0
    hcsl_sum = 0
    pcsl_sum = 0
    for m in np.arange(n):
        i = coords[0][m]
        j = coords[1][m]
        h_sum += hpsd[i,j]
        p_sum += ppsd[i,j]
        hcsl_sum += hcsl_ps[i,j]
        pcsl_sum += pcsl_ps[i,j]
    hpsd_fct[frq] = h_sum/n
    ppsd_fct[frq] = p_sum/n
    hcsl_psd_fct[frq] = hcsl_sum/n
    pcsl_psd_fct[frq] = pcsl_sum/n

neg = min(np.where(hcov_fct<0)[0])
Vh = hcov[0,0]
y = np.log(hcov_fct[1:neg]/Vh)
xih = -np.mean(d[1:neg]/y)      # fit exponential cov fct

neg = min(np.where(pcov_fct<0)[0])
Vp = hcov[0,0]
y = np.log(pcov_fct[1:neg]/Vp)
xip = -np.mean(d[1:neg]/y)      # fit exponential cov fct

neg = min(np.where(hcsl_cov_fct<0)[0])
Vh_csl = hcsl_cov[0,0]
y = np.log(hcsl_cov_fct[1:neg]/Vh_csl)
xih_csl = -np.mean(d[1:neg]/y)      # fit exponential cov fct
len(np.where(xihs>xih_csl)[0])/len(xihs) # interesting...

neg = min(np.where(pcsl_cov_fct<0)[0])
Vp_csl = pcsl_cov[0,0]
y = np.log(pcsl_cov_fct[1:neg]/Vp_csl)
xip_csl = -np.mean(d[1:neg]/y)      # fit exponential cov fct
len(np.where(xips>xip_csl)[0])/len(xips) # bummer

# visual inspection of cov fcts suggest nu = 0.5
# that is, cov fcts are decaying exp fcts
# and should be fit up to first negative value as tail end is garbage
plt.figure()
plt.scatter(x=d,y=hcov_fct)
plt.savefig('hcov_fct.png')
plt.close()
plt.figure()
plt.scatter(x=d,y=pcov_fct)
plt.savefig('pcov_fct.png')
plt.close()
plt.figure()
plt.scatter(x=d,y=hcsl_cov_fct)
plt.savefig('hcsl_cov_fct.png')
plt.close()
plt.figure()
plt.scatter(x=d,y=pcsl_cov_fct)
plt.savefig('pcsl_cov_fct.png')
plt.close()

plt.figure()
plt.scatter(x=d,y=hpsd_fct)
plt.savefig('hpsd_fct.png')
plt.close()
plt.figure()
plt.scatter(x=d,y=ppsd_fct)
plt.savefig('ppsd_fct.png')
plt.close()
plt.figure()
plt.scatter(x=d,y=hcsl_psd_fct)
plt.savefig('hcsl_psd_fct.png')
plt.close()
plt.figure()
plt.scatter(x=d,y=pcsl_psd_fct)
plt.savefig('pcsl_psd_fct.png')
plt.close()

# histograms of marg vars
plt.figure()
plt.hist(Vhs, bins=30, alpha=0.5, density=True, stacked=True)
plt.scatter(x=Vh_csl,y=0,c="red")
plt.savefig('Vh_hist.png')
plt.close()
plt.figure()
plt.hist(Vps, bins=30, alpha=0.5, density=True, stacked=True)
plt.scatter(x=Vp_csl,y=0,c="red")
plt.savefig('Vp_hist.png')
plt.close()

# histograms of char scales
plt.figure()
plt.hist(xihs, bins=30, alpha=0.5, density=True, stacked=True)
plt.scatter(x=xih_csl,y=0,c="red")
plt.savefig('xih_hist.png')
plt.close()
plt.figure()
plt.hist(xips, bins=30, alpha=0.5, density=True, stacked=True)
plt.scatter(x=xip_csl,y=0,c="red")
plt.savefig('xip_hist.png')
plt.close()

#
# now trying stuff w cross-psd and cross-cov... and w coh...
#

xps = np.zeros((len(h_normies)*len(p_normies),res,res//2+1), dtype=complex)
xcovs = np.zeros((len(h_normies)*len(p_normies),res,res))
i = 0
k = 0
for l_h in h_normies:
    j = 0
    for l_p in p_normies:
        xps[k,:,:], xcovs[k,:,:] = cf.fft_cov(h.p[i,:,:],p.p[j,:,:])
        j += 1
        k += 1
    i += 1
xpsd = np.mean(xps,axis=0)
xcov = np.fft.irfft2(xpsd)

csl_xps, csl_xcov = cf.fft_cov(h.p[h_normies[20],:,:],p.p[p_normies[10],:,:])

xcov_fct = np.zeros(len(d))
xpsd_fct = np.zeros(len(f))
csl_xcov_fct = np.zeros(len(d))
csl_xpsd_fct = np.zeros(len(f))
for dst in np.arange(len(d)):
    coords = np.where(dists==d[dst])
    n = len(coords[0])
    x_sum = 0
    csl_sum = 0
    for m in np.arange(n):
        i = coords[0][m]
        j = coords[1][m]
        x_sum += xcov[i,j]
        csl_sum += csl_xcov[i,j]
    xcov_fct[dst] = x_sum/n
    csl_xcov_fct[dst] = csl_sum/n
for frq in np.arange(len(f)): # replace with a build_fct function
    coords = np.where(freqs==f[frq])
    n = len(coords[0])
    x_sum = 0
    csl_sum = 0
    for m in np.arange(n):
        i = coords[0][m]
        j = coords[1][m]
        x_sum += abs(xpsd[i,j])
        csl_sum += abs(csl_xps[i,j])
    xpsd_fct[frq] = x_sum/n
    csl_xpsd_fct[frq] = csl_sum/n

# think this needs to be done before averaging...
csl_coh_fct = csl_xpsd_fct/np.sqrt(hpsd_fct*ppsd_fct)

plt.figure()
plt.scatter(x=d,y=csl_coh_fct)
plt.savefig('csl_coh_fct.png')
plt.close()

plt.figure()
plt.scatter(x=d,y=csl_xpsd_fct)
plt.savefig('csl_xpsd_fct.png')
plt.close()

plt.figure()
plt.scatter(x=d,y=xpsd_fct)
plt.savefig('xpsd_fct.png')
plt.close()

plt.figure()
plt.scatter(x=d,y=xcov_fct)
plt.savefig('xcov_fct.png')
plt.close()

plt.figure()
plt.scatter(x=d,y=csl_xcov_fct)
plt.savefig('csl_xcov_fct.png')
plt.close()

plt.figure()
matplt = plt.matshow(hcov)
plt.colorbar(matplt)
plt.savefig('h_cov.png')
plt.close()

plt.figure()
matplt = plt.matshow(hpsd)
plt.colorbar(matplt)
plt.savefig('hpsd.png')
plt.close()

plt.figure()
matplt = plt.matshow(pcov)
plt.colorbar(matplt)
plt.savefig('p_cov.png')
plt.close()

plt.figure()
matplt = plt.matshow(ppsd)
plt.colorbar(matplt)
plt.savefig('ppsd.png')
plt.close()

plt.figure()
matplt = plt.matshow(xcov)
plt.colorbar(matplt)
plt.savefig('x_cov.png')
plt.close()

plt.figure()
matplt = plt.matshow(csl_xcov)
plt.colorbar(matplt)
plt.savefig('csl_x_cov.png')
plt.close()

# repeat iLD analysis for different dispersal distances in each spp
#   - for each pair of dispersal dists, estimate spatial scale of allele freq turnover
#   - fit neutral iLD to beta distr for each spp separately
#   - estimate p-value of causal locus for each spp from resp beta distr
#   - plot p-values as function of spatial scale of allele freq turnover separately for ea spp
#   - repeat this whole thing for each pair of weak, moderate and strong biotic selection

#
# estimate selection coefficients
#

# first define model parameters (pulled from tree seq file)
ts = pyslim.load("mam.trees")
pars = cf.Pars(
    SI = ts.metadata["SLiM"]["user_metadata"]["SI"][0],         # spatial interspp intxn dist
    hc = ts.metadata["SLiM"]["user_metadata"]["hc"][0],       # multiplicative cost on host fitness for each infection
    pc = ts.metadata["SLiM"]["user_metadata"]["pc"][0],       # multiplicative cost on parasite fitness for failing to infect
    pb = ts.metadata["SLiM"]["user_metadata"]["pb"][0],       # multiplicative benefit on parasite fitness for successful infection
    minpr = ts.metadata["SLiM"]["user_metadata"]["minpr"][0],    # minimum probability of infection
    maxpr = ts.metadata["SLiM"]["user_metadata"]["maxpr"][0],    # maximum probability of infection
    lmbda = ts.metadata["SLiM"]["user_metadata"]["lmbda"][0],    # width of spatial niche
    ch = ts.metadata["SLiM"]["user_metadata"]["ch"][0],     # strength of host competition
    cp = ts.metadata["SLiM"]["user_metadata"]["cp"][0])     # strength of para competition

# compute global selection coefficients
h_s, p_s = cf.glbl_selCoeffs(sys, pars, h_pr=1.0, p_pr=1.0)
