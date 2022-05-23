import msprime, pyslim
import os

def add_tick_cycle(ts):
    tables = ts.dump_tables()
    ts_metadata = tables.metadata
    ts_metadata['SLiM']['cycle'] = 1
    ts_metadata['SLiM']['tick'] = 1
    tables.metadata = ts_metadata
    return tables.tree_sequence()

hts = msprime.sim_ancestry(samples=10, population_size=10, sequence_length=1e8)
hts = pyslim.annotate_defaults(hts, model_type="nonWF", slim_generation=1)
hts = add_tick_cycle(hts)

hts.dump(os.path.expanduser('~/gsccs-data/hinit.trees'))

pts = msprime.sim_ancestry(samples=10, population_size=10, sequence_length=1e8)
pts = pyslim.annotate_defaults(hts, model_type="nonWF", slim_generation=1)
pts = add_tick_cycle(hts)

pts.dump(os.path.expanduser('~/gsccs-data/pinit.trees'))

# here pop_0 will be an unsampled population that doesn't affect anything
demog = msprime.Demography.from_species_tree(
        "(pop_0:100000000.0,pop_1:100000000.0)A",
        initial_size={"pop_0": 10, "pop_1": 10, "A": 10}
)

pts = msprime.sim_ancestry(samples={"pop_0": 0, "pop_1": 10}, demography=demog, sequence_length=1e8)
# assert np.all(pts.tables.nodes.population == 1)
pts = pyslim.annotate_defaults(pts, model_type="nonWF", slim_generation=1)
pts = add_tick_cycle(pts)
