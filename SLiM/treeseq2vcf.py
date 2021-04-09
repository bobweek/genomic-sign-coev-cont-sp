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

#
# output vcf's
#

# positions and vcf files for the host (p1)
indivlist = []
indivnames = []
with open("host", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"]
                               + ["time", "x", "y", "z"]) + "\n") #this is whatever metadata you're interested in
    for i in range(p1_ts.num_individuals):
        indivlist.append(i)
        ind = p1_ts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(ind.time),
                str(ind.location[0]), str(ind.location[1]), str(ind.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("host.vcf", "w") as vcffile:
  p1_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)

# positions and vcf files for the parasite (p2)
indivlist = []
indivnames = []
with open("parasite", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"]
                               + ["time", "x", "y", "z"]) + "\n") #this is whatever metadata you're interested in
    for i in range(p2_ts.num_individuals):
        indivlist.append(i)
        ind = p2_ts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(ind.time),
                str(ind.location[0]), str(ind.location[1]), str(ind.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("parasite.vcf", "w") as vcffile:
  p2_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)

  # todo: compute spatial autocorrelation and cross-correlation functions in R!
