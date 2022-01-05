from os import chdir
from os import listdir
import matplotlib.pyplot as plt
import pyslim
import tskit
import msprime
import numpy as np
import pandas as pd
import pprint

#
# load in data
#

# L_dat = pd.read_csv('pars.txt', delimiter = ",") # needs work
L1 = 500000 # num neutral loci shandwhiching caus spp1, eventually read these in from file
L2 = 500000 # num neutral loci shandwhiching caus spp2

print(f"\n")
print(f"LOADING COMBINED TREE SEQ:")
slim_ts = pyslim.load("mam.trees")
print(f"The tree sequence has {slim_ts.num_trees} trees on a genome of length {slim_ts.sequence_length},"
      f" {slim_ts.num_individuals} individuals, {slim_ts.num_samples} 'sample' genomes,"
      f" and {slim_ts.num_mutations} mutations.")
print(f"\n")

# relable derived states to fit binaryMutation() model
tbls = slim_ts.dump_tables()
ads = []
for muts in slim_ts.tables.mutations:
      if len(muts.derived_state) > 0:
            ads.append('1')
      if len(muts.derived_state) == 0:
            ads.append('0')
tbls.mutations.packset_derived_state(ads)
aas = []
for syt in slim_ts.tables.sites:
      aas.append('0')
tbls.sites.packset_ancestral_state(aas)
# tbls.sort() # i think its already sorted
slim_ts = pyslim.SlimTreeSequence(tbls.tree_sequence())

#
# split the combined tree seq into spp specific tree seq's
#

print(f"SPLITTING TREE SEQ:")

# subset by species
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

# remove genetic information in adjacent species region
rm_intvl = np.array([2*L1+1,2*L1+2*L2+2])
rm_intvl = rm_intvl[np.newaxis,:]
h_ts = pyslim.SlimTreeSequence(h_ts.delete_intervals(rm_intvl))
h_ts = pyslim.SlimTreeSequence(h_ts.rtrim())
rm_intvl = np.array([0,2*L1+1])
rm_intvl = rm_intvl[np.newaxis,:]
p_ts = pyslim.SlimTreeSequence(p_ts.delete_intervals(rm_intvl))
p_ts = pyslim.SlimTreeSequence(p_ts.ltrim())

# print some stuff
print(f"The host tree sequence has {h_ts.num_trees} trees on a genome of length {h_ts.sequence_length},"
      f" {h_ts.num_individuals} individuals, {h_ts.num_samples} 'sample' genomes,"
      f" and {h_ts.num_mutations} mutations.")
print(f"The parasite tree sequence has {p_ts.num_trees} trees on a genome of length {p_ts.sequence_length},"
      f" {p_ts.num_individuals} individuals, {p_ts.num_samples} 'sample' genomes,"
      f" and {p_ts.num_mutations} mutations.")
print(f"\n")

#
# recapitation
#

print(f"RECAPITATING EACH SPP SEPARATELY:")

# perhaps set Ne to harmonic mean of historical pop sizes?
h_recap_ts = h_ts.recapitate(recombination_rate=1e-7, Ne=28000, random_seed=4)
p_recap_ts = p_ts.recapitate(recombination_rate=1e-7, Ne=36000, random_seed=5)

# sprinkle on mutations
h_ts = pyslim.SlimTreeSequence(
    msprime.sim_mutations(h_recap_ts, rate=1e-8, model=msprime.BinaryMutationModel()))
p_ts = pyslim.SlimTreeSequence(
    msprime.sim_mutations(p_recap_ts, rate=1e-8, model=msprime.BinaryMutationModel()))

# save recapped trees
h_ts.dump("h.mam.recap.trees")
p_ts.dump("p.mam.recap.trees")

# print some output
print(f"Host tree sequence now has {h_ts.num_trees} trees,"
      f" and {h_ts.num_mutations} mutations.")
print(f"Parasite tree sequence now has {p_ts.num_trees} trees,"
      f" and {p_ts.num_mutations} mutations.")
print(f"\n")

#
# output vcf's and metadata
#

print(f"SAVING VCFS AND METADATA:")

# finding the causal locus of the host
for i in h_ts.sites():
      if(i.position==L1):
            h_causL = i.id
# finding the causal locus of the parasite
for i in p_ts.sites():
      if(i.position==L1):
            p_causL = i.id

# write out site id's of causal loci for ea spp
with open("causL.mam.txt", "w") as causLfile:
      causLfile.writelines("\t".join(["h", "p"]) + "\n")
      causLfile.writelines("\t".join([str(h_causL), str(p_causL)]) + "\n") # (x,y)=coords, z=trait

# positions and vcf files for the host (h)
indivlist = []
indivnames = []
with open("host.mam.txt", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y", "z"]) + "\n") # (x,y)=coords, z=trait
    for i in h_ts.individuals():
        indivlist.append(i.id)
        vcf_label = f"tsk_{i.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(i.location[0]), str(i.location[1]), str(i.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("host.mam.vcf", "w") as vcffile:
      h_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)

# positions and vcf files for the parasite (p)
indivlist = []
indivnames = []
with open("parasite.mam.txt", "w") as indfile:
    indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y", "z"]) + "\n") # (x,y)=coords, z=trait
    for i in p_ts.individuals():
        indivlist.append(i.id)
        vcf_label = f"tsk_{i.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(i.location[0]), str(i.location[1]), str(i.location[2])]
        indfile.writelines("\t".join(data) + "\n")
with open("parasite.mam.vcf", "w") as vcffile:
      p_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)
