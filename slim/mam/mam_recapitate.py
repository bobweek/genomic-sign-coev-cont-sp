# import tskit
import pyslim
import msprime
import numpy as np

#
# load in data
#

print(f"\n")
print(f"Loading combined tree seq:")
ts = pyslim.load("mam.trees")
print(f"The tree sequence has {ts.num_trees} trees on a genome of length {ts.sequence_length},"
      f" {ts.num_individuals} individuals, {ts.num_samples} 'sample' genomes,"
      f" and {ts.num_mutations} mutations.")
print(f"\n")

# pull number of neutral loci surrounding causal locus for each spp
Lh = ts.metadata["SLiM"]["user_metadata"]["Lh"][0]
Lp = ts.metadata["SLiM"]["user_metadata"]["Lp"][0]

# relabel derived states to fit binaryMutation() model
tbls = ts.dump_tables()
ads = []
for muts in ts.tables.mutations:
      if len(muts.derived_state) > 0:
            ads.append('1')
      if len(muts.derived_state) == 0:
            ads.append('0')
tbls.mutations.packset_derived_state(ads)
aas = []
for syt in ts.tables.sites:
      aas.append('0')
tbls.sites.packset_ancestral_state(aas)
# tbls.sort() # i think its already sorted
ts = pyslim.SlimTreeSequence(tbls.tree_sequence())

# save memory
# del ts

#
# split the combined tree seq into spp specific tree seq's
#

print(f"Splitting tree seq by species:")

# subset by species
h_inds = np.where(ts.individual_populations == 1)[0]
h_nodes = []
for i in h_inds:
   h_nodes.extend(ts.individual(i).nodes)
h_ts = ts.simplify(h_nodes,keep_input_roots=True)
h_tbls = h_ts.tables

p_inds = np.where(ts.individual_populations == 2)[0]
p_nodes = []
for i in p_inds:
   p_nodes.extend(ts.individual(i).nodes)
p_ts = ts.simplify(p_nodes,keep_input_roots=True)
p_tbls = p_ts.tables

# subset by spp using tbls (is it worth it?)
# h_tbls = tskit.TableCollection(sequence_length=2*Lh+1)
# for n in tbls.nodes:
#       if n.population == 1 and n.time == 0.0:
#             h_tbls.nodes.add_row(n) // don't work

# remove genetic information in adjacent species region
rm_intvl = np.array([2*Lh+1,2*Lh+2*Lp+2])
rm_intvl = rm_intvl[np.newaxis,:]
h_tbls.delete_intervals(rm_intvl)
h_tbls.rtrim()
h_ts = pyslim.SlimTreeSequence(h_tbls.tree_sequence())

rm_intvl = np.array([0,2*Lh+1])
rm_intvl = rm_intvl[np.newaxis,:]
p_tbls.delete_intervals(rm_intvl)
p_tbls.ltrim()
p_ts = pyslim.SlimTreeSequence(p_tbls.tree_sequence())

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

print(f"Recapitating each species separately:")

# perhaps set Ne to harmonic mean of historical pop sizes?
h_Ne = h_ts.num_individuals
p_Ne = p_ts.num_individuals
h_recap_ts = h_ts.recapitate(recombination_rate=1e-7, Ne=h_Ne, random_seed=4)
p_recap_ts = p_ts.recapitate(recombination_rate=1e-7, Ne=p_Ne, random_seed=5)

# sprinkle on mutations
h_ts = pyslim.SlimTreeSequence(
    msprime.sim_mutations(h_recap_ts, rate=1e-8, model=msprime.BinaryMutationModel()))
p_ts = pyslim.SlimTreeSequence(
    msprime.sim_mutations(p_recap_ts, rate=1e-8, model=msprime.BinaryMutationModel()))

# save recapped trees
# h_ts.dump("h.mam.recap.trees")
# p_ts.dump("p.mam.recap.trees")

# print some output
print(f"Host tree sequence now has {h_ts.num_trees} trees,"
      f" and {h_ts.num_mutations} mutations.")
print(f"Parasite tree sequence now has {p_ts.num_trees} trees,"
      f" and {p_ts.num_mutations} mutations.")
print(f"\n")

#
# output GenotypeArrays and metadata
#

print(f"Saving GenotypeArrays (in scikit-allel format) and metadata:")

# saving snp locations and
# finding the causal locus of the host
h_snps = []
for i in h_ts.sites():
      h_snps.append(i.position)
      if(i.position==Lh):
            h_causL = i.id
np.savetxt("h_snps.csv", h_snps, delimiter=",")

# saving snp locations and
# finding the causal locus of the parasite
p_snps = []
for i in p_ts.sites():
      p_snps.append(i.position)
      if(i.position==Lh):
            p_causL = i.id
np.savetxt("p_snps.csv", p_snps, delimiter=",")

# write out site id's of causal loci for ea spp
with open("causL.mam.txt", "w") as causLfile:
      causLfile.writelines("\t".join(["h", "p"]) + "\n")
      causLfile.writelines("\t".join([str(h_causL), str(p_causL)]) + "\n") # (x,y)=coords, z=trait

# positions for the host (h)
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
gm = h_ts.genotype_matrix()
gm1 = gm[:,0:h_ts.num_individuals]
gm2 = gm[:,h_ts.num_individuals:h_ts.num_samples]
ga = np.reshape(np.array([gm1,gm2]),(h_ts.num_sites,h_ts.num_individuals,2))
np.save("host.mam",ga)

# positions for the parasite (p)
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
gm = p_ts.genotype_matrix()
gm1 = gm[:,0:p_ts.num_individuals]
gm2 = gm[:,p_ts.num_individuals:p_ts.num_samples]
ga = np.reshape(np.array([gm1,gm2]),(p_ts.num_sites,p_ts.num_individuals,2))
np.save("para.mam",ga)