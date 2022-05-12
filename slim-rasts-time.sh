#!/bin/bash

slim slim/multispp-slim/para-host.slim
./rasts.sh
./time-series.sh

ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'slimulation rasters and analysis complete!'