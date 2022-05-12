#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate slim-tskit
rm ~/gsccs-data/rast-data/*.csv
python python/rasts.py
rm ~/gsccs-data/rast-data/*.png
R -f r/rasts.R
cd ~/gsccs-data/rast-data/
rm *.mp4
ffmpeg -r 20 -i N%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p N.mp4 -y
ffmpeg -r 20 -i z%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p z.mp4 -y
ffmpeg -r 20 -i v%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p v.mp4 -y
ffmpeg -r 20 -i Nz%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p Nz.mp4 -y

ntfy -b pushover -o user_key uagcx5q4jqpv21t1mowsp5dfecjut8 send 'rasters complete!'