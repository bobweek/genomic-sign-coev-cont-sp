import matplotlib.pyplot as plt
import pyslim
import tskit
import msprime
import numpy as np
import pandas
import pprint

# is it possible to run slim from python????

L1 = 10  # num qtl for each spp
L2 = 1e4 # num neutral loci

slim_ts = pyslim.load(
      "/home/bb/gits/genomic-sign-coev-cont-sp/slim/coev_sim.trees")
print(f"The tree sequence has {slim_ts.num_trees} trees on a genome of length {slim_ts.sequence_length},"
      f" {slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,"
      f" and {slim_ts.num_mutations} mutations.")

# split the tree sequence
# need to also subset genome by qtls
# right now i think the spare qtls in each seq are treated as neutral loci?
h_inds = np.where(slim_ts.individual_populations == 1)[0]
h_nodes = []
for i in h_inds:
   h_nodes.extend(slim_ts.individual(i).nodes)
h_ts = slim_ts.simplify(h_nodes,keep_input_roots=True)
p_inds = np.where(slim_ts.individual_populations == 2)[0]
p_nodes = []
for i in p_inds:
   p_nodes.extend(slim_ts.individual(i).nodes)
p_ts = slim_ts.simplify(p_nodes,keep_input_roots=True)

# h_locs = slim_ts.individual_locations[h_inds]
# p_locs = slim_ts.individual_locations[p_inds]

# yay! they're the same!!!
h_ts.num_nodes + p_ts.num_nodes == slim_ts.num_nodes
h_ts.num_individuals + p_ts.num_individuals == slim_ts.num_individuals

# old approach probably doesn't work right
# slim_ts.individuals_alive_at(0)
# slim_ts.num_individuals
# p_inds = np.where(slim_ts.individual_populations == 2)[0]
# h_nodes = []
# for i in h_inds:
#       h_nodes.extend(slim_ts.individual(i).nodes)
# h_ts = pyslim.SlimTreeSequence(slim_ts.subset(h_nodes))
# p_nodes = []
# for i in p_inds:
#       p_nodes.extend(slim_ts.individual(i).nodes)
# p_ts = pyslim.SlimTreeSequence(slim_ts.subset(p_nodes))

# perform recapitation
# gotta for new pyslim to use recomb map
# recomb_pos = [0, 2*L1-1, 2*L1+L2-1]
# recomb_rates = [0.5, 1e-5, 0]
# recomb_map = msprime.RecombinationMap(positions = recomb_pos, rates = recomb_rates)
# h_recap_ts = h_ts.recapitate(recombination_map=recomb_map, Ne=10000) # perhaps set Ne to harmonic mean of historical pop sizes?
# p_recap_ts = p_ts.recapitate(recombination_map=recomb_map, Ne=10000)

h_recap_ts = h_ts.recapitate(recombination_rate=1e-8, Ne=10000) # perhaps set Ne to harmonic mean of historical pop sizes?
p_recap_ts = p_ts.recapitate(recombination_rate=1e-8, Ne=10000)

# sprinkle on mutations
# is this using the correct mutation model??
h_ts = pyslim.SlimTreeSequence(
    msprime.mutate(h_recap_ts, rate=1e-5, keep=True))
p_ts = pyslim.SlimTreeSequence(
    msprime.mutate(p_recap_ts, rate=1e-5, keep=True))

# save recapped trees
h_ts.dump("/home/bb/gits/genomic-sign-coev-cont-sp/slim/h.recap.trees")
p_ts.dump("/home/bb/gits/genomic-sign-coev-cont-sp/slim/p.recap.trees")

# print some output
print(f"Species 1 tree sequence now has {h_ts.num_trees} trees,"
      f" and {h_ts.num_mutations} mutations.")
print(f"Species 2 tree sequence now has {p_ts.num_trees} trees,"
      f" and {p_ts.num_mutations} mutations.")

#
# plotting locations and coloring by trait values
#

# locations of species one alive at end of simulation
h_locs = h_ts.individual_locations

# locations of species two alive at end of simulation
p_locs = p_ts.individual_locations

# plot spatial locations and color by trait for both species
fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("Species One")
ax.scatter(h_locs[:, 0], h_locs[:, 1], s=10, c=h_locs[:, 2])
ax = fig.add_subplot(122)
ax.set_title("Species Two")
ax.scatter(p_locs[:, 0], p_locs[:, 1], s=10, c=p_locs[:, 2])
fig.savefig("/home/bb/gits/genomic-sign-coev-cont-sp/slim/coev_sim_locations.png")

#
# output vcf's for most recent time (t=0)
#

# positions and vcf files for the host (h)
indivlist = []
indivnames = []
with open("host.txt", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y", "z"]) + "\n") #this is whatever metadata you're interested in
    for i in h_ts.individuals:
        indivlist.append(i)
        ind = h_ts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(ind.location[0]), str(ind.location[1]), str(ind.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("host.vcf", "w") as vcffile:
      h_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)

# positions and vcf files for the parasite (p)
indivlist = []
indivnames = []
with open("parasite.txt", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y", "z"]) + "\n") #this is whatever metadata you're interested in
    for i in p_ts.individuals:
        indivlist.append(i)
        ind = p_ts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(ind.location[0]), str(ind.location[1]), str(ind.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("parasite.vcf", "w") as vcffile:
      p_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)

# union treeseqs to output combined vcf
# also output txt file with ind ID, spp, location and trait value
both_ts = pyslim.SlimTreeSequence(h_ts.union(p_ts, np.repeat(tskit.NULL,p_ts.num_nodes), add_populations=True))

indivlist = []
indivnames = []
with open("both.txt", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"]
                               + ["species", "x", "y", "z"]) + "\n") #this is whatever metadata you're interested in
    for i in both_ts.individuals_alive_at(0):
        indivlist.append(i)
        ind = both_ts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(ind.population), str(ind.location[0]), str(ind.location[1]), str(ind.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("both.vcf", "w") as vcffile:
      both_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)
