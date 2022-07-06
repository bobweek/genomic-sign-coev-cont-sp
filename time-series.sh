#!/bin/bash

Rscript r/make-ts-report.R
Rscript -e "rsconnect::deployApp('r/time-series')" -y
# ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'time-series report deployed!'