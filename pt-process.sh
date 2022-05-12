#!/bin/bash

rm ~/gsccs-data/ind-data/*.png
R -f r/pt-plts.R
ffmpeg -r 30 -i ~/gsccs-data/ind-data/z%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p r/pt-process.mp4 -y
ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'animation complete!'