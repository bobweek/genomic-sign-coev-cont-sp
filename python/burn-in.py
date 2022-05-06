from statistics import mean
import msprime
import pyslim
import numpy as np

N = 10000
n = 100
L = 1e8
mu = 1e-13
k = 1
rho = 1e-8

# first simulate whole genome to get num segr sites
# to get freely recomb loci, then run a sim per segr site

# generate genealogy
# demog_model = msprime.Demography() # this was producing two populations per species
# demog_model.add_population(initial_size=N)
# hts = msprime.sim_ancestry(samples=n,demography=demog_model,sequence_length=L,recombination_rate=rho)
# pts = msprime.sim_ancestry(samples=n,demography=demog_model,sequence_length=L,recombination_rate=rho)
hts = msprime.sim_ancestry(samples=n,population_size=N,sequence_length=L,recombination_rate=rho)
pts = msprime.sim_ancestry(samples=n,population_size=N,sequence_length=L,recombination_rate=rho)

# annotate
hts = pyslim.annotate_defaults(hts, model_type="nonWF", slim_generation=1)
pts = pyslim.annotate_defaults(pts, model_type="nonWF", slim_generation=1)

# slim mutations
hmut = msprime.SLiMMutationModel(type=1) # type = 1 for m1 (spp 1)
pmut = msprime.SLiMMutationModel(type=2) # type = 1 for m1 (spp 1)
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

# check we didn't mess anything up
assert htables.mutations.num_rows == hts.num_mutations
assert ptables.mutations.num_rows == pts.num_mutations
print(f"The host trait effects range from {min(hmut.values()):0.2e}"f" to {max(hmut.values()):0.2e}.")
print(f"The para trait effects range from {min(pmut.values()):0.2e}"f" to {max(pmut.values()):0.2e}.")

# what kind slim we runnin here??
htables.metadata
ptables.metadata

ptables.populations[0]

# gotta change a few things
ts_metadata = htables.metadata
ts_metadata["SLiM"]["spatial_dimensionality"] = "xy"
htables.metadata = ts_metadata
htables.populations.clear()
for p in hts.populations():
    pm = p.metadata
    pm['slim_id'] = 0
    pm['bounds_x1'] = 100
    pm['bounds_y1'] = 100
    htables.populations.add_row(metadata=pm)
hts = htables.tree_sequence()
hts.dump("/home/bb/gsccs-data/hinit.trees")
ts_metadata = ptables.metadata
ts_metadata["SLiM"]["spatial_dimensionality"] = "xy"
ptables.metadata = ts_metadata
ptables.populations.clear()
ptables.populations.add_row()
for p in pts.populations(): # make empty pop in pop[0]
    pm = p.metadata
    pm['slim_id'] = 1
    pm['bounds_x1'] = 100
    pm['bounds_y1'] = 100
    ptables.populations.add_row(metadata=pm)
pts = ptables.tree_sequence()
pts.dump("/home/bb/gsccs-data/pinit.trees")

# trait values
Zh = np.zeros(n)
Zp = np.zeros(n)
hgm = hts.genotype_matrix()
pgm = pts.genotype_matrix()
for i in 2*np.arange(n):
    Zh[int(i/2)] = sum(hgm[:,i]*hfx)+sum(hgm[:,i+1]*hfx)
    Zp[int(i/2)] = sum(pgm[:,i]*pfx)+sum(pgm[:,i+1]*pfx)

np.var(Zh,ddof=1)
np.var(Zp,ddof=1)
