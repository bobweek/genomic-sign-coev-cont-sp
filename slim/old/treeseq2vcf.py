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


# todo: compute spatial autocorrelation and cross-correlation functions in R!

#
# combine recapitated tree seqs in one big mommi tree seq
# this part is no longer needed, but may be useful later. will eventually place in different file 
#

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
