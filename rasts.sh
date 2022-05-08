#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate slim-tskit
python python/rasts.py
R -f r/rasts.R
cd ~/gsccs-data/rasts-data/
ffmpeg -r 20 -i N%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p N.mp4
ffmpeg -r 20 -i z%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p z.mp4
ffmpeg -r 20 -i v%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p v.mp4
ffmpeg -r 20 -i Nz%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p Nz.mp4