#!/bin/bash

rm ~/gsccs-data/ind-data/*.csv
slim slim/para-host.slim
ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'slimulation complete!'
./time-series.sh