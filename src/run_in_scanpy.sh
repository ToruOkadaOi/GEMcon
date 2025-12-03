#!/bin/bash

if command -v micromamba &> /dev/null; then # command -v checks if micromamba is installed
    micromamba run -n scanpy_legacy python "$@"
else
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate scanpy_legacy
    python "$@"
fi