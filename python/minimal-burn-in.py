from statistics import mean
import msprime
import pyslim
import numpy as np
import os

N = 1000
n = 100
L = 1e8
mu = 1e-13
k = 1 # size of gaussian mutations
rho = 1e-8

hts = msprime.sim_ancestry(samples=n,population_size=N,sequence_length=L,recombination_rate=rho)
pts = msprime.sim_ancestry(samples=n,population_size=N,sequence_length=L,recombination_rate=rho)

hts = pyslim.annotate_defaults(hts, model_type="nonWF", slim_generation=1)
pts = pyslim.annotate_defaults(pts, model_type="nonWF", slim_generation=1)

# slim mutations
hmut = msprime.SLiMMutationModel(type=1) # type = 1 for m1 (spp 1)
pmut = msprime.SLiMMutationModel(type=2) # type = 2 for m2 (spp 2)
hts = msprime.sim_mutations(hts,rate=mu,model=hmut,keep=True)
pts = msprime.sim_mutations(pts,rate=mu,model=pmut,keep=True)
print(f"The host now has {hts.num_mutations} mutations, at "f"{hts.num_sites} distinct sites.")
print(f"The para now has {pts.num_mutations} mutations, at "f"{pts.num_sites} distinct sites.")

# adding selection coefficients (which are trait effects for me)
htables = hts.tables
htables.mutations.clear()
hmut = {}
hfx = np.zeros(hts.num_sites)
for m in hts.mutations():
  md_list = m.metadata["mutation_list"]
  slim_ids = m.derived_state.split(",")
  assert len(slim_ids) == len(md_list)
  for sid, md in zip(slim_ids, md_list):
     if sid not in hmut:
        hfx[int(sid)] = np.random.normal(loc=0,scale=k) # gaussian mutations
        hmut[sid] = hfx[int(sid)]
     md["selection_coeff"] = hmut[sid]
  _ = htables.mutations.append(
          m.replace(metadata={"mutation_list": md_list})
  )
ptables = pts.tables
ptables.mutations.clear()
pmut = {}
pfx = np.zeros(pts.num_sites)
for m in pts.mutations():
  md_list = m.metadata["mutation_list"]
  slim_ids = m.derived_state.split(",")
  assert len(slim_ids) == len(md_list)
  for sid, md in zip(slim_ids, md_list):
     if sid not in pmut:
        pfx[int(sid)] = np.random.normal(loc=0,scale=k) # gaussian mutations
        pmut[sid] = pfx[int(sid)]
     md["selection_coeff"] = pmut[sid]
  _ = ptables.mutations.append(
          m.replace(metadata={"mutation_list": md_list})
  )

#
# editing spp1 metadata
#

ts_metadata = htables.metadata
ts_metadata["SLiM"]["spatial_dimensionality"] = "xy"
ts_metadata['SLiM']['cycle'] = 1
ts_metadata['SLiM']['tick'] = 1
htables.metadata = ts_metadata

# hindividual_metadata = [ind.metadata for ind in htables.individuals]
# for md in hindividual_metadata:
#    md["subpopulation"] = 0
#    ims = htables.individuals.metadata_schema
#    htables.individuals.packset_metadata(
#       [ims.validate_and_encode_row(md) for md in hindividual_metadata])

# hmut_metadata = [mut.metadata for mut in htables.mutations]
# for md in hmut_metadata:
#    md["mutation_list"][0]["subpopulation"] = 0
#    ims = htables.mutations.metadata_schema
#    htables.mutations.packset_metadata(
#       [ims.validate_and_encode_row(md) for md in hmut_metadata])

htables.populations.clear()
for p in hts.populations():
    pm = p.metadata
    pm['slim_id'] = 0
    pm['bounds_x1'] = 100
    pm['bounds_y1'] = 100
    htables.populations.add_row(metadata=pm)
hts = htables.tree_sequence()
hts.dump(os.path.expanduser('~/gsccs-data/hinit.trees'))

#
# editing spp1 metadata
#

ts_metadata = ptables.metadata
ts_metadata["SLiM"]["spatial_dimensionality"] = "xy"
# ts_metadata['SLiM']['cycle'] = 1
# ts_metadata['SLiM']['tick'] = 1
ptables.metadata = ts_metadata

# pindividual_metadata = [ind.metadata for ind in ptables.individuals]
# for md in pindividual_metadata:
#    md["subpopulation"] = 1
#    ims = ptables.individuals.metadata_schema
#    ptables.individuals.packset_metadata(
#       [ims.validate_and_encode_row(md) for md in pindividual_metadata])

# pmut_metadata = [mut.metadata for mut in ptables.mutations]
# for md in pmut_metadata:
#    md["mutation_list"][0]["subpopulation"] = 1
#    ims = ptables.mutations.metadata_schema
#    ptables.mutations.packset_metadata(
#       [ims.validate_and_encode_row(md) for md in pmut_metadata])

ptables.populations.clear()
ptables.populations.add_row()
for p in pts.populations(): # make empty pop in pop[0]
    pm = p.metadata
    pm['slim_id'] = 1
    pm['bounds_x1'] = 100
    pm['bounds_y1'] = 100
    ptables.populations.add_row(metadata=pm)
pts = ptables.tree_sequence()
pts.dump(os.path.expanduser('~/gsccs-data/pinit.trees'))
