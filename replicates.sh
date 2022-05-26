#!/bin/bash

cd peqg22
python replicates.py
ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'replicates finished!'