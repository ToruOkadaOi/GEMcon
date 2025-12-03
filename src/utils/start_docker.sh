#!/bin/bash

echo 'Insert directory to be binded to container(complete path):'

read d_path

docker run -d -p 8783:8787 -p 8899:8899 \
  -v $d_path:/home/rstudio \
  toluene123/scrna_complete

echo 'r-studio will be on localhost:8783 and jupyterlab on localhost:8899'

container=$(docker ps | grep 'toluene123/scrna_complete' | awk '{print $1}')

docker exec -it $container \
  jupyter notebook --ip=0.0.0.0 --port=8899 --no-browser --allow-root --NotebookApp.token='' --NotebookApp.password=''

echo done
