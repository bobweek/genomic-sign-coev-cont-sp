#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate slim-tskit
python python/slimulation-analysis.py
R -f r/visualize-time-series.R

ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'slimulation analysis complete!'