import matplotlib.pyplot as plt
import pyslim
import tskit
import msprime
import numpy as np
import pandas
import pprint

# load in slim tree seqs
p1_ts = pyslim.load(
      "/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/p1.recap.trees")
p2_ts = pyslim.load(
      "/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/p2.recap.trees")

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

# p1 & p2 alive at t=0 referenced in unioned tree seq
pp1_t0 = np.intersect1d(p_ts.individuals_alive_at(0), p1_ts.individuals_alive_at(0))
pp2_t0 = np.intersect1d(p_ts.individuals_alive_at(0), len(
    p1_ts.individuals_alive_at(0))+p2_ts.individuals_alive_at(0))

# sanity check
p_ts.pairwise_diversity(pp1_t0) == p1_ts.pairwise_diversity(
    p1_ts.individuals_alive_at(0))
p_ts.pairwise_diversity(pp2_t0) == p2_ts.pairwise_diversity(
    p2_ts.individuals_alive_at(0)) # why not working?


#
# trying to relabel populations by species id...
#

pprint.pprint(p1_ts.metadata_schema)

tflag = p1_ts.tables.nodes.flags
ttime = p1_ts.tables.nodes.time
tind  = p1_ts.tables.nodes.individual
p1_ts.tables.metadata_schema

p1_ts.tables.nodes[0]

encoded_metadata_column = [
    p1_tbls.metadata_schema.validate_and_encode_row(r) for r in metadata_column
]
metadata, metadata_offset = tskit.pack_bytes(encoded_metadata_column)

tmeta = p1_ts.tables.nodes.metadata

P1 = np.repeat(1, p1_ts.num_nodes)
p1_tbls.nodes.set_columns(
    flags=tflag, time=ttime, population=P1.astype(np.int32), individual=tind, metadata=tmeta)

tflag = p2_ts.tables.nodes.flags
ttime = p2_ts.tables.nodes.time
tind  = p2_ts.tables.nodes.individual
tmeta = p2_ts.tables.nodes.metadata
p2_tbls.nodes.set_columns(
    flags=tflag, time=ttime, population=np.repeat(1, p2_ts.num_nodes), individual=tind, metadata=tmeta)

#
# plotting locations and coloring by trait values
#

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
fig.savefig("/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/coev_sim_locations.png")

# todo: compute spatial autocorrelation and cross-correlation functions
