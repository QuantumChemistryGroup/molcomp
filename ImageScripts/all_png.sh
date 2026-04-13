#!/bin/sh

mogrify -trim $*  # Crop images
montage $* -geometry 300x+10+10 -tile 5x -pointsize 16 -set label '%t' all.png
pngquant -f -o all.png all.png # Compress image
rm $*
