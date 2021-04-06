import pyslim
import tskit
import msprime
import numpy as np
import pandas

slim_ts = pyslim.load(
      "/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/coev_sim.trees")
print(f"The tree sequence has {slim_ts.num_trees} trees on a genome of length {slim_ts.sequence_length},"
      f" {slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,"
      f" and {slim_ts.num_mutations} mutations.")

# split the tree sequence
p1_nodes = []
p2_nodes = []
for i in range(slim_ts.num_nodes):
      if slim_ts.tables.nodes[i].population == 1:
            p1_nodes.append(i)
      if slim_ts.tables.nodes[i].population == 2:
            p2_nodes.append(i)
p1_ts = pyslim.SlimTreeSequence(slim_ts.subset(p1_nodes))   # for species 1
p2_ts = pyslim.SlimTreeSequence(slim_ts.subset(p2_nodes))   # for species 2

# perform recapitation
p1_recap_ts = p1_ts.recapitate(recombination_rate=1e-8, Ne=10000)
p2_recap_ts = p2_ts.recapitate(recombination_rate=1e-8, Ne=10000)

# sprinkle on mutations
p1_ts = pyslim.SlimTreeSequence(
    msprime.mutate(p1_recap_ts, rate=1e-8, keep=True))
p2_ts = pyslim.SlimTreeSequence(
    msprime.mutate(p2_recap_ts, rate=1e-8, keep=True))

# make sure each tree has a single root
sum([t.num_roots == 1 for t in p1_ts.trees()]) == len(p1_ts.trees())
sum([t.num_roots == 1 for t in p2_ts.trees()]) == len(p2_ts.trees())

# save recapped trees
p1_ts.dump("/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/p1.recap.trees")
p2_ts.dump("/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/p2.recap.trees")

# print some output
print(f"Species 1 tree sequence now has {p1_ts.num_trees} trees,"
      f" and {p1_ts.num_mutations} mutations.")
print(f"Species 2 tree sequence now has {p2_ts.num_trees} trees,"
      f" and {p2_ts.num_mutations} mutations.")