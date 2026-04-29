#!/bin/bash
if command -v micromamba &> /dev/null; then # command -v checks if micromamba is installed
    # On Docker
    micromamba run -n gecko_aman python "$@"
else
    # without Docker
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate gecko_aman
    python "$@"
fi

### Setup gecko py conda env