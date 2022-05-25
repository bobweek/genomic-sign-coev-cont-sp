#
# this script loads in tree seqs outputted by slim
# exports individual locations and trait values to csv's
# then exports genotype arrays of causal mutations
# then wipes the causal mutations from the tree seq
# then sprinkles on neutral mutations
# then exports genotype arrays of neutral mutations
#

import os
import tskit
import msprime
import numpy as np


# load trees

host_ts = tskit.load(os.path.expanduser('~/gsccs-data/host-slim.trees'))
para_ts = tskit.load(os.path.expanduser('~/gsccs-data/para-slim.trees'))


# save causal snp locations

h_csl_snps = []
for i in host_ts.sites():
      h_csl_snps.append(i.position)
np.savetxt(os.path.expanduser('~/gsccs-data/h-csl-snps.csv'), h_csl_snps, delimiter=",")

p_csl_snps = []
for i in para_ts.sites():
      p_csl_snps.append(i.position)
np.savetxt(os.path.expanduser('~/gsccs-data/p-csl-snps.csv'), p_csl_snps, delimiter=",")


# export locations and trait values

indivlist = []
indivnames = []
with open(os.path.expanduser('~/gsccs-data/host.txt'), "w") as indfile:
    indfile.writelines(",".join(["vcf_label"]
                               + ["x", "y", "z"]) + "\n") # (x,y)=coords, z=trait
    for i in host_ts.individuals():
        indivlist.append(i.id)
        vcf_label = f"tsk_{i.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(i.location[0]), str(i.location[1]), str(i.location[2])]
        indfile.writelines(",".join(data) + "\n")

indivlist = []
indivnames = []
with open(os.path.expanduser('~/gsccs-data/para.txt'), "w") as indfile:
    indfile.writelines(",".join(["vcf_label"]
                               + ["x", "y", "z"]) + "\n") # (x,y)=coords, z=trait
    for i in para_ts.individuals():
        indivlist.append(i.id)
        vcf_label = f"tsk_{i.id}"
        indivnames.append(vcf_label)
        data = [vcf_label, str(i.location[0]), str(i.location[1]), str(i.location[2])]
        indfile.writelines(",".join(data) + "\n")


# export individual allele freqs and genotype matrices of causal mutations

hgm = host_ts.genotype_matrix()
hgm1 = hgm[:,0:host_ts.num_individuals]
hgm2 = hgm[:,host_ts.num_individuals:host_ts.num_samples]
hfrq = (hgm1+hgm2)/2
np.save(os.path.expanduser('~/gsccs-data/hfrq-causal'),hfrq)
hga = np.reshape(np.array([hgm1,hgm2]),(host_ts.num_sites,host_ts.num_individuals,2))
np.save(os.path.expanduser('~/gsccs-data/hga-causal'),hga)

pgm = para_ts.genotype_matrix()
pgm1 = pgm[:,0:para_ts.num_individuals]
pgm2 = pgm[:,para_ts.num_individuals:para_ts.num_samples]
pfrq = (pgm1+pgm2)/2
np.save(os.path.expanduser('~/gsccs-data/pfrq-causal'),pfrq)
pga = np.reshape(np.array([pgm1,pgm2]),(para_ts.num_sites,para_ts.num_individuals,2))
np.save(os.path.expanduser('~/gsccs-data/pga-causal'),pga)


# inspect trees

print(f"The host tree sequence has {host_ts.num_trees} trees on a genome of length {host_ts.sequence_length},"
      f" {host_ts.num_individuals} individuals, {host_ts.num_samples} 'sample' genomes,"
      f" and {host_ts.num_mutations} causal mutations.")

print(f"The tree sequence has {para_ts.num_trees} trees on a genome of length {para_ts.sequence_length},"
      f" {para_ts.num_individuals} individuals, {para_ts.num_samples} 'sample' genomes,"
      f" and {para_ts.num_mutations} causal mutations.")


# sprinkle on neutral mutations

mu = 1e-11 # neutral mutations occur an order of magnitude faster than causal

htables = host_ts.tables
htables.mutations.clear()
hts = htables.tree_sequence()

ptables = para_ts.tables
ptables.mutations.clear()
pts = ptables.tree_sequence()

hmut = msprime.SLiMMutationModel(type=1)
hts = msprime.sim_mutations(host_ts,rate=mu,model=hmut,keep=True)

pmut = msprime.SLiMMutationModel(type=1)
pts = msprime.sim_mutations(pts,rate=mu,model=pmut,keep=True)


# how many mutations now?

print(f"Host tree sequence has {hts.num_mutations} neutral mutations.")

print(f"Para tree sequence has {pts.num_mutations} neutral mutations.")


# save neutral snp locations

h_ntl_snps = []
for i in hts.sites():
      h_ntl_snps.append(i.position)
np.savetxt(os.path.expanduser('~/gsccs-data/h-ntl-snps.csv'), h_ntl_snps, delimiter=",")

p_ntl_snps = []
for i in pts.sites():
      p_ntl_snps.append(i.position)
np.savetxt(os.path.expanduser('~/gsccs-data/p-ntl-snps.csv'), p_ntl_snps, delimiter=",")

# export individual allele freqs and genotype matrices of neutral mutations

hgm = hts.genotype_matrix()
hgm1 = hgm[:,0:hts.num_individuals]
hgm2 = hgm[:,hts.num_individuals:hts.num_samples]
hfrq = (hgm1+hgm2)/2
np.save(os.path.expanduser('~/gsccs-data/hfrq-neutrl'),hfrq)
hga = np.reshape(np.array([hgm1,hgm2]),(hts.num_sites,hts.num_individuals,2))
np.save(os.path.expanduser('~/gsccs-data/hga-neutrl'),hga)

pgm = pts.genotype_matrix()
pgm1 = pgm[:,0:pts.num_individuals]
pgm2 = pgm[:,pts.num_individuals:pts.num_samples]
pfrq = (pgm1+pgm2)/2
np.save(os.path.expanduser('~/gsccs-data/pfrq-neutrl'),pfrq)
pga = np.reshape(np.array([pgm1,pgm2]),(pts.num_sites,pts.num_individuals,2))
np.save(os.path.expanduser('~/gsccs-data/pga-neutrl'),pga)