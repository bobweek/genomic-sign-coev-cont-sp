import pyslim
import tskit
import msprime
import numpy as np
import pandas

p1_ts = pyslim.load(
      "/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/p1.recap.trees")
p2_ts = pyslim.load(
      "/home/bb/Projects/The Genomic Signature of Coevolution in Continuous Space/SLiM/p2.recap.trees")

p1_ts.coiterate(p2_ts)
