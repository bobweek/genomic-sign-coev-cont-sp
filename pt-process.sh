#!/bin/bash

R -f r/pt-plts.R
ffmpeg -r 20 -i ~/gsccs-data/ind-data/z%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ~/gsccs-data/pt-process.mp4 -y
cd ~/software
./gdrive upload ~/gsccs-data/pt-process.mp4
# for f in ~/gsccs-data/ind-data/*.png; do ./gdrive upload $f; done
# ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'animation complete!'