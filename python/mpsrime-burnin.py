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
