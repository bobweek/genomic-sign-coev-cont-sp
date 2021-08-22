import matplotlib.pyplot as plt
from numpy.core.records import array
import pyslim
import tskit
import msprime
import numpy as np
import pandas as pd
import pprint

space_width = 100/2

h_ts = pyslim.load(
      "/home/bb/gits/genomic-sign-coev-cont-sp/slim/h.recap.trees")
p_ts = pyslim.load(
      "/home/bb/gits/genomic-sign-coev-cont-sp/slim/p.recap.trees")

both_ts = pyslim.SlimTreeSequence(h_ts.union(p_ts, np.repeat(tskit.NULL,p_ts.num_nodes), add_populations=True))

# bin space into K x K grid

K = 20
binwidth = 2*space_width/K
bins = range()

h_ts.individual_locations[:,range(2)].max()

# draw smaller samples of inds for plotting
sSize = 100
ht0 = np.random.choice(range(h_ts.num_individuals),sSize,replace=False)
pt0 = np.random.choice(range(p_ts.num_individuals),sSize,replace=False)

hn0 = []
hid = []
for i in ht0:
    hn0.append(h_ts.individual(i).nodes)
    hid.append(i)
pn0 = []
pid = []
for i in pt0:
    pn0.append(p_ts.individual(i).nodes)
    pid.append(i)

h_pairs = [(i, j) for i in range(sSize) for j in range(i, sSize)]
h_div = h_ts.divergence(hn0, indexes=h_pairs)
h_div_sep = h_ts.divergence(hn0, windows=[0,10,20,h_ts.sequence_length], indexes=h_pairs)
h_q_div = h_div_sep[0]
h_n_div = h_div_sep[2]


p_pairs = [(i, j) for i in range(sSize) for j in range(i, sSize)]
p_div = p_ts.divergence(pn0, indexes=p_pairs)
p_div_sep = p_ts.divergence(pn0, windows=[0,10,20,p_ts.sequence_length], indexes=p_pairs)
p_q_div = p_div_sep[1]
p_n_div = p_div_sep[2]

#
# collect arrays
# correcting for periodic boundaries so that max geo dist = sqrt(2)*space_width
#

h_geog_dist = np.repeat(0.0, len(h_pairs))
h_phen_dist = np.repeat(0.0, len(h_pairs))
h_locs = h_ts.individual_locations[:,range(2)][ht0]
h_phen = h_ts.individual_locations[:,2][ht0]
for k, (i, j) in enumerate(h_pairs):
   h_geog_dist[k] = np.sqrt(np.sum((abs(h_locs[i] - h_locs[j])%space_width)**2))
   h_phen_dist[k] = np.abs(h_phen[i] - h_phen[j])

# h_phen_pd = pd.Series(h_ts.individual_locations[:,2])
# h_phen_pd.plot.hist(grid=True, bins=10, rwidth=0.9,
#                    color='#607c8e')
# plt.title('Distribution Host Trait')
# plt.xlabel('Trait Value')
# plt.ylabel('Count')
# plt.grid(axis='y', alpha=0.5)
# plt.show()

p_geog_dist = np.repeat(0.0, len(p_pairs))
p_phen_dist = np.repeat(0.0, len(p_pairs))
p_locs = p_ts.individual_locations[:,range(2)]
p_phen = p_ts.individual_locations[:,2]
for k, (i, j) in enumerate(p_pairs):
   p_geog_dist[k] = np.sqrt(np.sum((abs(p_locs[pid[i], :] - p_locs[pid[j], :])%space_width)**2))
   p_phen_dist[k] = np.abs(p_phen[pid[i]] - p_phen[pid[j]])

# genetic dist X geo dist
fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("Host Species")
ax.scatter(h_geog_dist, 1e3 * h_n_div, s=20, alpha=0.5)
ax.set_xlabel("geographic distance")
ax.set_ylabel("genetic distance (diffs/Kb)")
ax = fig.add_subplot(122)
ax.set_title("Parasite Species")
ax.scatter(p_geog_dist, 1e3 * p_n_div, s=20, alpha=0.5)
ax.set_xlabel("geographic distance")
ax.set_ylabel("genetic distance (diffs/Kb)")
fig.savefig("/home/bb/gits/genomic-sign-coev-cont-sp/slim/spatial_sim_ibd.png")

# genetic dist X geo dist
# using qtls only
fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("Host Species")
ax.scatter(h_geog_dist, h_q_div, s=20, alpha=0.5)
ax.set_xlabel("geographic distance")
ax.set_ylabel("genetic distance (just qtls)")
ax = fig.add_subplot(122)
ax.set_title("Parasite Species")
ax.scatter(p_geog_dist, p_q_div, s=20, alpha=0.5)
ax.set_xlabel("geographic distance")
ax.set_ylabel("genetic distance (just qtls)")
fig.savefig("/home/bb/gits/genomic-sign-coev-cont-sp/slim/spatial_sim_ibd_q.png")

# genetic dist X pheno dist
# using qtls only
fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("Host Species")
ax.scatter(h_phen_dist, h_q_div, s=20, alpha=0.5)
ax.set_xlabel("phenotypic distance")
ax.set_ylabel("genetic distance (diffs/Kb)")
ax = fig.add_subplot(122)
ax.set_title("Parasite Species")
ax.scatter(p_phen_dist, p_q_div, s=20, alpha=0.5)
ax.set_xlabel("phenotypic distance")
ax.set_ylabel("genetic distance (diffs/Kb)")
fig.savefig("/home/bb/gits/genomic-sign-coev-cont-sp/slim/spatial_sim_ibp.png")

# pheno dist X geo dist
fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("Host Species")
ax.scatter(h_geog_dist, h_phen_dist, s=20, alpha=0.5)
ax.set_xlabel("geographic distance")
ax.set_ylabel("phenotypic distance")
ax = fig.add_subplot(122)
ax.set_title("Parasite Species")
ax.scatter(p_geog_dist, p_phen_dist, s=20, alpha=0.5)
ax.set_xlabel("geographic distance")
ax.set_ylabel("phenotypic distance")
fig.savefig("/home/bb/gits/genomic-sign-coev-cont-sp/slim/spatial_sim_pbd.png")

#
# todo w/in tskit: 
#
# compute phenotypic auto/cross-correlation functions and compare to spde model <- THIS!
#

# first divide space into equally sized cells (say 100 units to begin with)

div = 10
cell_width = 2*space_width/div
cell_locs = []
cell_h_means = []
cell_p_means = []
for i in range(div):
    for j in range(div):
        x = (i+0.5)*cell_width
        y = (j+0.5)*cell_width
        cell_locs.append([x,y])
        hx = h_ts.individual_locations[:,0]
        hy = h_ts.individual_locations[:,1]
        px = p_ts.individual_locations[:,0]
        py = p_ts.individual_locations[:,1]
        h_within = np.logical_and(np.logical_and(hx <= x + 0.5*cell_width, hx >= x - 0.5*cell_width), np.logical_and(hy <= y + 0.5*cell_width, hy >= y - 0.5*cell_width))
        p_within = np.logical_and(np.logical_and(px <= x + 0.5*cell_width, px >= x - 0.5*cell_width), np.logical_and(py <= y + 0.5*cell_width, py >= y - 0.5*cell_width))
        hz = np.mean(h_ts.individual_locations[h_within,2])
        pz = np.mean(p_ts.individual_locations[p_within,2])
        cell_h_means.append(hz)
        cell_p_means.append(pz)

# look at correlation in discretized space (how does this depend on scale?)
# looks like the pattern is not very striking, which motivates the need for genomic method
# plt.scatter(cell_h_means, cell_p_means, s=20, alpha=0.5)
# plt.show()

hZ = np.mean(cell_h_means)
pZ = np.mean(cell_p_means)

# next step is to compute correlation as fct of distance

# first compute an array comprising all possible distances (modulo space_width)
# then loop through cells and for each cell look at all other cells of equal distance away and average product of difference to mean
#   average these averages for each cell iterated through


# how do allele freqs correlate with trait in other species?
# how do allele freqs correlate with trait of same species?