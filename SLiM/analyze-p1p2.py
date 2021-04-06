import matplotlib.pyplot as plt
import pyslim
import tskit
import msprime
import numpy as np
import pandas

# load in slim tree seqs
p1_ts = pyslim.load(
      "/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/p1.recap.trees")
p2_ts = pyslim.load(
      "/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/p2.recap.trees")

# get tables for each spp
p1_tbls = p1_ts.dump_tables()
p2_tbls = p2_ts.dump_tables()

# make unioned table
p_tbls = p1_ts.dump_tables()
p_tbls.union(p2_tbls, np.repeat(-1, p2_ts.num_nodes))

# convert to unioned slim tree seq
p_ts = pyslim.SlimTreeSequence(p_tbls.tree_sequence())

# sanity check
p_ts.num_nodes == p1_ts.num_nodes + p2_ts.num_nodes

# population ids for unioned tree seq 
# this doesn't work since unioned tree lists alive at t=0 first, regardless of population
#p1_ids = np.arange(p1_ts.num_nodes)
#p2_ids = p1_ts.num_nodes + np.arange(p2_ts.num_nodes)

# analysis
p1_ts.pairwise_diversity(p1_ts.individuals_alive_at(0))
p2_ts.pairwise_diversity(p2_ts.individuals_alive_at(0))
p_ts.pairwise_diversity(p_ts.individuals_alive_at(0))

# p1 & p2 alive at t=0 referenced in unioned tree seq
pp1_t0 = np.intersect1d(p_ts.individuals_alive_at(0), p1_ts.individuals_alive_at(0))
pp2_t0 = np.intersect1d(p_ts.individuals_alive_at(0), len(
    p1_ts.individuals_alive_at(0))+p2_ts.individuals_alive_at(0))

# sanity check
p_ts.pairwise_diversity(pp1_t0) == p1_ts.pairwise_diversity(
    p1_ts.individuals_alive_at(0))
p_ts.pairwise_diversity(pp2_t0) == p2_ts.pairwise_diversity(
    p2_ts.individuals_alive_at(0)) # why not working?

p1_ts.Fst(p1_ts.individuals_alive_at(0))

# subset of species one alive at end of simulation
p1_t0 = p1_ts.individuals_alive_at(0)
# locations of species one alive at end of simulation
p1_locs = p1_ts.individual_locations[p1_t0,:]

# subset of species two alive at end of simulation
p2_t0 = p2_ts.individuals_alive_at(0)
# locations of species two alive at end of simulation
p2_locs = p2_ts.individual_locations[p2_t0, :]

# plot spatial locations and color by trait for both species
fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("Species One")
ax.scatter(p1_locs[:, 0], p1_locs[:, 1], s=10, c=p1_locs[:, 2])
ax = fig.add_subplot(122)
ax.set_title("Species Two")
ax.scatter(p2_locs[:, 0], p2_locs[:, 1], s=10, c=p2_locs[:, 2])
fig.savefig("/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/coev_sim_locations.png")

p1_ts.tables.nodes[1]

p2_ts.tables.nodes[1].population = 2

# todo: compute spatial autocorrelation and cross-correlation functions
