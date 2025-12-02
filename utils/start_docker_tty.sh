#!/bin/bash

read d_path
d_path=${d_path:-$(pwd)}

docker run --rm -it -p 8783:8787 -p 8899:8899 \
  -v "$d_path":/home/rstudio \
  toluene123/scrna_complete /bin/bash
