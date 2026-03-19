#!/bin/bash
## assumes you have a cplex lic.; install accordingly to the env named cplex
cp -r /opt/ibm/ILOG/CPLEX_Studio2211/cplex/python/3.10/x86-64_linux /tmp/cplex_python
micromamba run -n cplex_aman_new pip install /tmp/cplex_python
rm -rf /tmp/cplex_python
