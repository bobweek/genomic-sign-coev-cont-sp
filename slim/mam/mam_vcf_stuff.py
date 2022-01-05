import allel
import numpy as np
import skgstat as skg
import matplotlib.pyplot as plt 
import mam_classes_fcts as cf
# from sklearn.neighbors import KernelDensity
# from sklearn.model_selection import GridSearchCV
from scipy.stats import t, normaltest, kstest, anderson, shapiro, binned_statistic

# for reloading custom module while developing
import importlib
# importlib.reload(cf) # example use

#
# load in data
#

# yields type(sys) = System
sys = cf.loadUp('host.mam.txt','parasite.mam.txt','causL.mam.txt','host.mam.vcf','parasite.mam.vcf')

#
# drop low frequency loci (prevents large spurious correlations when we discretize)
#

# yields type(sys) = System
sys = cf.dropSmolFreqs(sys, p_min=0.2)

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

#
# discretize space into res x res grid of cells
#

# spatial resolution
res = 32

# spatial extent
width = 100 # eventually read in from file
height = 100 # eventually read in from file

# discretization, yields type(dis_sys) = discSystem
dis_sys = cf.discSpace(sys, width, height, res)
h, p = dis_sys

# prevent unnecessary use of memory
# del sys

#
# plot allele freq surfaces of causal loci and three randomly selected neutral loci
#

cf.plot_p(dis_sys,res)

# #
# # filter loci using single spp Fst patterns (hold off for now)
# #

# how does host-para coev affect Fst?
h_Fst, p_Fst = cf.Fst(dis_sys)

# make "manhattan" plots
cf.Fst_plts(h_Fst,p_Fst,dis_sys)

# coev appears to significantly lower host Fst
# coev appears to slightly increase para Fst
# this seems to agree with the surface plots
# does this pattern generalize?
# for what region of parameter space?
# what other patterns occur?

# filter out host loci in highest 90% of Fst scores
dis_sys.h = cf.singleSppFilter(-h_Fst, 90, h)

# filter out para loci in lowest 90% of Fst scores
dis_sys.p = cf.singleSppFilter(p_Fst, 90, p)

# update species specific objects
h, p = dis_sys

#
# compute interspecific spatial correlations of allele frequencies
#

# compute
iscaf = cf.discaf(dis_sys)

# fit a t-distr
t_stats = []
for l_h in np.arange(h.S):
    for l_p in np.arange(p.S):
        r = iscaf[l_h,l_p]
        t_stt = r*np.sqrt((2.5*res-2)/(1-r**2))
        t_stats.append(t_stt)

# see if computed t-statistics fit the student t-distr
# using kolmogorov-smirnov test
kstest(t_stats,"t",args=(2*res-2,))

# visual check if t is good fit
plt.figure()
xrang = np.arange(min(t_stats),max(t_stats),0.001)
plt.scatter(x=xrang,y=t.pdf(xrang,2.5*res-2),c="red")
plt.scatter(x=np.mean(t_stats),y=0,c="black")
plt.hist(t_stats, bins=100, alpha=0.5, density=True, stacked=True)
plt.savefig('t_hist_clr.png')
plt.close()

# are the causal loci correlated more than any other pair?
# cf.causal_search(iscaf, dis_sys)

# make manhattan plots for causal loci
cf.csl_mnhtn_plts(abs(iscaf),dis_sys)

#
# distribution of iscaf across independent (w/in spp) loci
#

# pull out iscaf across unlinked loci
iscaf_unlinked = abs(iscaf[h_unlinked,:][:,p_unlinked])

# make histogram of iscaf at unlinked loci and superimpose causal iscaf
# eventually use cross-val to find best bandwidth
# grid = GridSearchCV(KernelDensity(),
#                     {'bandwidth': np.linspace(0.1, 5.0, 30)},
#                     cv=20) # 20-fold cross-validation
# grid.fit(np.array(iscaf_unlinked).reshape(-1,1))
# print grid.best_params_
cf.iscaf_plts(iscaf,abs(iscaf[h.causL,p.causL]))

#
# does student's t-distr apply?
#

# check if allele frequencies are normally distributed across locations
# using three different tests
h_shaps = []
h_andrs = []
h_nrmlt = []
for l in h_unlinked:
    h_shaps.append(shapiro(h.p[l,:,:].flatten()).pvalue)
    h_andrs.append(anderson(h.p[l,:,:].flatten()).significance_level)
    h_nrmlt.append(normaltest(h.p[l,:,:].flatten()).pvalue)
p_shaps = []
p_andrs = []
p_nrmlt = []
for l in p_unlinked:
    p_shaps.append(shapiro(p.p[l,:,:].flatten()).pvalue)
    p_andrs.append(anderson(p.p[l,:,:].flatten()).significance_level)
    p_nrmlt.append(normaltest(p.p[l,:,:].flatten()).pvalue)

# just using shapiro-wilks test for now
h_normies = np.where(np.array(h_shaps) < 0.05)[0]
p_normies = np.where(np.array(p_shaps) < 0.05)[0]

# proportion of loci that are non-normal
# worth finding a transform to try and normalize (logit)?
1-len(h_normies)/len(h_shaps)
1-len(p_normies)/len(p_shaps)

#
# removing intraspecific spatial correlations of allele freqs
# so we can apply students t
#

# estimate covariance function
# this assumes allele freq surfs are sampled from same distribution
# so we can average spatial power spectra across loci to obtain cov
hpmean = np.mean(h.p.flatten())
hps = np.zeros((len(h_normies),res,res//2+1))
k = 0
for l in h_normies:
    hfft = np.fft.rfft2(h.p[l,:,:]-hpmean)
    hps[k,:,:] = abs(hfft)**2/res**2
    k += 1
hpsd = np.mean(hps,axis=0)
hcov = np.fft.irfft2(hpsd)

ppmean = np.mean(p.p.flatten())
pps = np.zeros((len(p_normies),res,res//2+1))
k = 0
for l in p_normies:
    pfft = np.fft.rfft2(p.p[l,:,:]-ppmean)
    pps[k,:,:] = abs(pfft)**2/res**2
    k += 1
ppsd = np.mean(pps,axis=0)
pcov = np.fft.irfft2(ppsd)

# once the covariance function is obtained need to compute expctd covariance
# of allele freqs at a single locus w/in a spp between every pair of cells
# yielding a res**2 x res**2 matrix of allele freq correlations
coords = []
for i in np.arange(res):
    for j in np.arange(res):
        coords.append([i,j])
hpcov = np.zeros((res**2,res**2))
ppcov = np.zeros((res**2,res**2))
coords = np.array(coords)
for i in np.arange(res**2):
    for j in np.arange(res**2):
        x = (coords[i]-coords[j])%res
        hpcov[i,j] = hcov[x[0],x[1]]
        ppcov[i,j] = pcov[x[0],x[1]]

# compute cholesky of inverse cov matrix
hpchol = np.linalg.cholesky(np.linalg.inv(hpcov))
ppchol = np.linalg.cholesky(np.linalg.inv(ppcov))

# decorrelate
hp_ucnr = np.zeros((len(h_normies),res**2))
k = 0
for l in h_normies:
    hp_ucnr[k] = np.matmul(np.transpose(hpchol),h.p[l,:,:].flatten())
    k += 1
pp_ucnr = np.zeros((len(p_normies),res**2))
k = 0
for l in p_normies:
    pp_ucnr[k] = np.matmul(np.transpose(ppchol),p.p[l,:,:].flatten())
    k += 1

# plot covariance function
# oh crap they seem to be highly correlated everywhere
plt.figure()
matplt = plt.matshow(hcov)
plt.colorbar(matplt)
plt.savefig('h_cov.png')
plt.close()
plt.figure()
matplt = plt.matshow(pcov)
plt.colorbar(matplt)
plt.savefig('p_cov.png')
plt.close()

# compute effective sample sizes

# marginal variances, aka nuggets
Vh = hcov[0,0]
Vp = pcov[0,0]

# variance of sample means across loci
hp_means = []
for l_h in h_normies:
    hp_means.append(np.mean(h.p[l_h,:,:]))
Vh_hat = np.var(np.array(hp_means))
pp_means = []
for l_p in p_normies:
        pp_means.append(np.mean(p.p[l_p,:,:]))
Vp_hat = np.var(np.array(pp_means))

# effective sample sizes
neff_h = Vh/Vh_hat
neff_p = Vp/Vp_hat
df = neff_h*neff_p-2

# for the loci that are normally distributed
# compute student t-stat
t_stats = []
for l_h in np.arange(len(h_normies)):
    for l_p in np.arange(len(p_normies)):
        hp = hp_ucnr[l_h]
        pp = pp_ucnr[l_p]
        r = np.corrcoef(hp,pp)[0,1]
        t_stt = r*np.sqrt((res**2-2)/(1-r**2))
        t_stats.append(t_stt)

# see if computed t-statistics fit the student t-distr
# using kolmogorov-smirnov test
kstest(t_stats,"t",args=(res**2-2,))

# visual check if t is good fit
plt.figure()
xrang = np.arange(min(t_stats),max(t_stats),0.001)
plt.scatter(x=xrang,y=t.pdf(xrang,res**2-2),c="red")
plt.scatter(x=np.mean(t_stats),y=0,c="black")
plt.hist(t_stats, bins=100, alpha=0.5, density=True, stacked=True)
plt.savefig('t_hist.png')
plt.close()

#
# when does spatial autocorrelation mask the iLD signature of coevolution?
#

# repeat iLD analysis for different dispersal distances in each spp
#   - for each pair of dispersal dists, estimate spatial scale of allele freq turnover
#   - fit neutral iLD to beta distr for each spp separately
#   - estimate p-value of causal locus for each spp from resp beta distr
#   - plot p-values as function of spatial scale of allele freq turnover separately for ea spp
#   - repeat this whole thing for each pair of weak, moderate and strong biotic selection

#
# how does spatial scale for causal locus compare to that of neutral loci?
# can we use this as a single species filter to reduce number of interspp comparisons?
# UNDER CONSTRUCTION
#

# compute distr of spatial scales across loci w/in ea spp separately
# is spatial scale of causal locus significantly:
#   - longer? shorter?
# does this depend on relative dispersal abilities?
# can this be used as a single spp filter? (Fst might not be so reliable)


#
# PCA (have no idea if this is useful)
#

h_coords, h_model = allel.pca(hgt_sqsh)

p_coords, p_model = allel.pca(pgt_sqsh)

plt.figure()
plt.scatter(h_coords[:,0],h_coords[:,1])
plt.savefig('h_pca.png')
plt.close()

plt.figure()
plt.scatter(p_coords[:,0],p_coords[:,1])
plt.savefig('p_pca.png')
plt.close()

#
# estimate selection coefficients
#

# first define model parameters (eventually pull from file)
pars = cf.Pars(
    SI=0.5,         # spatial interspp intxn dist
    hc = 0.8,       # multiplicative cost on host fitness for each infection
    pc = 0.6,       # multiplicative cost on parasite fitness for failing to infect
    pb = 1.3,       # multiplicative benefit on parasite fitness for successful infection
    minpr = 0.1,    # minimum probability of infection
    maxpr = 1.0,    # maximum probability of infection
    lmbda = 0.5,    # width of spatial niche
    ch = 0.005,     # strength of host competition
    cp = 0.005)     # strength of para competition

# compute global selection coefficients
h_s, p_s = cf.glbl_selCoeffs(sys, pars, h_pr=1.0, p_pr=1.0)
