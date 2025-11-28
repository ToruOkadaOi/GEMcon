#!/bin/bash
if command -v micromamba &> /dev/null; then # command -v checks if micromamba is installed
    # On Docker
    micromamba run -n cplex_aman python "$@"
else
    # without Docker
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate cplex_aman_new
    python "$@"
fi
