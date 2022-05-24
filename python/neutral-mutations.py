# import matplotlib.pyplot as plt
import os
import tskit
import msprime
import numpy as np

import pandas
import pprint

# neutral mutations occur an order of magnitude faster than causal
mu = 1e-11

# load trees

host_ts = tskit.load(os.path.expanduser('~/gsccs-data/host.trees'))
para_ts = tskit.load(os.path.expanduser('~/gsccs-data/para.trees'))

# export genotype matrices of causal mutations

pgm = para_ts.genotype_matrix()
pgm1 = pgm[:,0:para_ts.num_individuals]
pgm2 = pgm[:,para_ts.num_individuals:para_ts.num_samples]
pga = np.reshape(np.array([pgm1,pgm2]),(para_ts.num_sites,para_ts.num_individuals,2))
np.save(os.path.expanduser('~/gsccs-data/pga-causal'),pga)


# inspect trees

print(f"The host tree sequence has {host_ts.num_trees} trees on a genome of length {host_ts.sequence_length},"
      f" {host_ts.num_individuals} individuals, {host_ts.num_samples} 'sample' genomes,"
      f" and {host_ts.num_mutations} mutations.")

print(f"The tree sequence has {para_ts.num_trees} trees on a genome of length {para_ts.sequence_length},"
      f" {para_ts.num_individuals} individuals, {para_ts.num_samples} 'sample' genomes,"
      f" and {para_ts.num_mutations} mutations.")

# sprinkle on mutations

hmut = msprime.SLiMMutationModel(type=3, next_id=host_ts.num_mutations) # type = 3 for m3 (neutral type for spp 1)
hts = msprime.sim_mutations(host_ts,rate=mu,model=hmut,keep=True)

pmut = msprime.SLiMMutationModel(type=4, next_id=para_ts.num_mutations) # type = 4 for m4 (neutral type for spp 2)
pts = msprime.sim_mutations(para_ts,rate=mu,model=pmut,keep=True)

# how many mutations now?

print(f"Species 1 tree sequence now has {hts.num_trees} trees,"
      f" and {hts.num_mutations} mutations.")

print(f"Species 2 tree sequence now has {pts.num_trees} trees,"
      f" and {pts.num_mutations} mutations.")

# save mutated trees

hts.dump(os.path.expanduser('~/gsccs-data/h-mut.trees'))

pts.dump(os.path.expanduser('~/gsccs-data/p-mut.trees'))

#
# output vcf's
#

hts.genotype_matrix()

pgm = pts.genotype_matrix()
pgm1 = pgm[:,0:pts.num_individuals]
pgm2 = pgm[:,pts.num_individuals:pts.num_samples]
pga = np.reshape(np.array([pgm1,pgm2]),(pts.num_sites,pts.num_individuals,2))
np.save(os.path.expanduser('~/gsccs-data/pga'),pga)


with open(os.path.expanduser('~/gsccs-data/host.vcf'), "w") as vcffile:
    hts.write_vcf(vcffile)

with open(os.path.expanduser('~/gsccs-data/para.vcf'), "w") as vcffile:
    pts.write_vcf(vcffile)
