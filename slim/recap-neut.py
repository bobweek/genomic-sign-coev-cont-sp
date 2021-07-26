import matplotlib.pyplot as plt
import pyslim
import tskit
import msprime
import numpy as np
import pandas
import pprint

L = 1e5 # num neutral loci

slim_ts = pyslim.load(
      "/home/bb/gits/genomic-sign-coev-cont-sp/slim/neut_sim.trees")
print(f"The tree sequence has {slim_ts.num_trees} trees on a genome of length {slim_ts.sequence_length},"
      f" {slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,"
      f" and {slim_ts.num_mutations} mutations.")

recap_ts = slim_ts.recapitate(recombination_rate=1e-8, Ne=10000) # perhaps set Ne to harmonic mean of historical pop sizes?

alive_inds = recap_ts.individuals_alive_at(0)
simp_ts = recap_ts.simplify(alive_inds)

# sprinkle on mutations
# is this using the correct mutation model??
ts = pyslim.SlimTreeSequence(
    msprime.mutate(simp_ts, rate=1e-5, keep=True))

# save recapped trees
ts.dump("/home/bb/gits/genomic-sign-coev-cont-sp/slim/neut.recap.trees")

# print some output
print(f"Species 1 tree sequence now has {ts.num_trees} trees,"
      f" and {ts.num_mutations} mutations.")

#
# plotting locations and coloring by trait values
#

# locations of species one alive at end of simulation
locs = ts.individual_locations

# plot spatial locations and color by trait for both species
fig = plt.figure(figsize=(6, 6), dpi=300)
ax = fig.add_subplot(111)
ax.set_title("Location")
ax.scatter(locs[:, 0], locs[:, 1], s=10)
fig.savefig("/home/bb/gits/genomic-sign-coev-cont-sp/slim/neut_sim_locations.png")

#
# output vcf's for most recent time (t=0)
#

# positions and vcf files
indivlist = []
indivnames = []
with open("/home/bb/gits/genomic-sign-coev-cont-sp/slim/neut.txt", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"] + ["x", "y", "z"]) + "\n") #this is whatever metadata you're interested in
    for i in ts.individuals():
        indivlist.append(i.id)
        # ind = ts.individual(i)
        vcf_label = f"tsk_{i.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(i.location[0]), str(i.location[1]), str(i.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("/home/bb/gits/genomic-sign-coev-cont-sp/slim/neut.vcf", "w") as vcffile:
      ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)
