import allel
import os
import numpy as np
import pandas as pd
import statistics as st
import matplotlib.pyplot as plt 
from skgstat import Variogram

plt.style.use('ggplot')

space_width = 100/2

# read in locations and trait values for each sample
# species = 0 corresponds to host
# species = 1 corresponds to parasite
metadat = pd.read_csv('/home/bb/gits/genomic-sign-coev-cont-sp/slim/both.txt', delimiter = "\t")

# number of samples for each species
h_samps = sum(metadat['species']==0)
p_samps = sum(metadat['species']==1)

# include part where indices of causal and neutral loci are read-in
# host_causaL = read in from file
# host_neutrL = read in from file
# para_causaL = read in from file
# para_neutrL = read in from file

# read in vcf files
h_callset = allel.read_vcf(os.path.expanduser('~/gsccs-data/host.vcf'))
pre_h_gt = allel.GenotypeArray(h_callset['calldata/GT'])

p_callset = allel.read_vcf('/home/bb/gits/genomic-sign-coev-cont-sp/slim/parasite.vcf')
pre_p_gt = allel.GenotypeArray(p_callset['calldata/GT'])

# "number of segregating sites"
# needs number of causal loci subtracted off
pre_h_S = pre_h_gt.n_variants

pre_p_S = pre_p_gt.n_variants

#
# subtracting off each others causal loci
#

# host:
# causal loci for host = 0:9
# causal loci for parasite = 10:19
hLrang1 = [x for x in range(10)]
hLrang2 = [x for x in range(20,pre_h_S)]
hLrang = np.union1d(hLrang1,hLrang2)
h_tri_gt = pre_h_gt.subset(sel0=hLrang)
# h_gt = allel.GenotypeArray(h_tri_gt % 2) # don't do these as the 0,1,2 are needed for diploids
h_gt = allel.GenotypeArray(h_tri_gt)
h_alt = h_gt.to_n_alt()
h_ac = h_gt.count_alleles()

# parasite
pLrang = [x for x in range(10,pre_p_S)]
p_tri_gt = pre_p_gt.subset(sel0=pLrang)
# p_gt = allel.GenotypeArray(p_tri_gt % 2) # don't do these as the 0,1,2 are needed for diploids
p_gt = allel.GenotypeArray(p_tri_gt)
p_alt = p_gt.to_n_alt()
p_ac =  p_gt.count_alleles()

# actual number of segregating sites
h_S = pre_h_S - 10
p_S = pre_p_S - 10

# pull out locations and trait values
h_inds = metadat['species']==0
p_inds = metadat['species']==1
h_loc = metadat[['x','y']][h_inds]
p_loc = metadat[['x','y']][p_inds]
h_trt = metadat['z'][h_inds]
p_trt = metadat['z'][p_inds]

#
# calculating statistics
#

# variogram for trait values

# how to force smoothness parameter of matern function to 1.0?

hV = Variogram(coordinates=h_loc, values=h_trt, model='matern', dist_func=lambda u, v: np.sqrt(((np.abs(u-v)%space_width)**2).sum()))
pV = Variogram(coordinates=p_loc, values=p_trt, model='matern', dist_func=lambda u, v: np.sqrt(((np.abs(u-v)%space_width)**2).sum()))

print(hV)
print(pV)

hV.plot()
pV.plot()

hV.describe()['params']

# calculating pi for hosts, parasites and combined
h_pi = allel.sequence_diversity(range(h_S), h_ac)
p_pi = allel.sequence_diversity(range(p_S), p_ac)

#
# find patterns of intra/interspecific LD in the neutral*neutral and causal*causal comparisons
#

h_R = allel.rogers_huff_r(h_alt)
p_R = allel.rogers_huff_r(p_alt)
allel.plot_pairwise_ld(h_R)
# plt.show()

# what is this doing?
hp_R = allel.rogers_huff_r_between(h_alt,p_alt)

#
# todo w/in scikit-allel:
#
# look a
#