#!/bin/bash

python python/burn-in.py
slim slim/mortal-para-host.slim
ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'slimulation complete!'