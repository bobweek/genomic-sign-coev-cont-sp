#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate slim-tskit
python python/slimulation-analysis.py
R -f r/visualize-time-series.R