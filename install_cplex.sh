#!/bin/bash
# Install CPLEX into the conda environment
# You need your own IBM CPLEX license — https://www.ibm.com/products/ilog-cplex-optimization-studio

# Set CPLEX_PATH to your installation location
# Default: /opt/ibm/ILOG/CPLEX_Studio2211 (standard IBM install path)
# Override with: export CPLEX_PATH=/your/path/to/CPLEX_Studio2211
CPLEX_PATH=${CPLEX_PATH:-/opt/ibm/ILOG/CPLEX_Studio2211}

if [ ! -d "$CPLEX_PATH/cplex" ]; then
    echo "CPLEX not found at $CPLEX_PATH"
    echo "Please ensure CPLEX is installed and set CPLEX_PATH if needed"
    echo "  export CPLEX_PATH=/path/to/your/CPLEX_Studio2211"
    exit 1
fi

echo "Installing CPLEX from $CPLEX_PATH ..."

cp -r "$CPLEX_PATH/cplex/python/3.10/x86-64_linux" /tmp/cplex_python
micromamba run -n cplex_aman_new pip install /tmp/cplex_python
rm -rf /tmp/cplex_python

echo "CPLEX installed successfully!"
