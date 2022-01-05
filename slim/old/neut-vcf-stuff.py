import allel
import numpy as np
import pandas as pd
import statistics as st
import matplotlib.pyplot as plt 
from rpy2.robjects.packages import importr

geoR = importr('geoR')

# from skgstat import Variogram

plt.style.use('ggplot')

space_width = 1/2

# read in locations and trait values for each sample
# species = 0 corresponds to host
# species = 1 corresponds to parasite
metadat = pd.read_csv('/home/bb/gits/genomic-sign-coev-cont-sp/slim/neut.txt', delimiter = "\t")

# number of samples for each species
# samps = sum(metadat['species']==0)

# include part where indices of causal and neutral loci are read-in
# host_causaL = read in from file
# host_neutrL = read in from file
# para_causaL = read in from file
# para_neutrL = read in from file

# read in vcf files
callset = allel.read_vcf('/home/bb/gits/genomic-sign-coev-cont-sp/slim/neut.vcf')
pre_gt = allel.GenotypeArray(callset['calldata/GT'])

gt = allel.GenotypeArray(pre_gt)
alt = gt.to_n_alt()
ac =  gt.count_alleles()

# "number of segregating sites"
S = pre_gt.n_variants

# pull out locations and trait values
# inds = range(0,metadat["x"].size)
loc = metadat[['x','y']][0:metadat["x"].size]
trt = metadat['z'][0:metadat["x"].size]

#
# calculating statistics
#

# variogram for trait values

# how to force smoothness parameter of matern function to 1.0?

V = Variogram(coordinates=loc, values=trt, model='matern', dist_func=lambda u, v: np.sqrt(((np.abs(u-v)%space_width)**2).sum()))

print(V)
V.plot()
V.describe()['params']

# calculating pi for hosts, parasites and combined
pi = allel.sequence_diversity(range(S), ac)

#
# find patterns of intra/interspecific LD in the neutral*neutral and causal*causal comparisons
#

R = allel.rogers_huff_r(alt)
allel.plot_pairwise_ld(R)
plt.show()

#
# todo w/in scikit-allel:
#
# look a
#