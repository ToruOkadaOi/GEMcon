# GEMcon

[![test_ga](https://github.com/ToruOkadaOi/GEMcon/actions/workflows/test.yml/badge.svg)](https://github.com/ToruOkadaOi/GEMcon/actions/workflows/test.yml) [![Docker](https://img.shields.io/docker/v/toluene123/gemcon?label=toluene123/gemcon)](https://hub.docker.com/r/toluene123/gemcon)

Pipeline for fetching expression /proteomic abundance data and building **context-specific human metabolic models** ([Human-GEM](https://github.com/SysBioChalmers/Human-GEM) based)

---

## Overview

GEMcon integrates data from:
- scRNA-seq ([Human Cell Atlas](https://www.humancellatlas.org))
- bulk RNA-seq ([GTEx](https://gtexportal.org/home/)) -->proposed
- proteomics ([PaxDB](https://pax-db.org))

---

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/ToruOkadaOi/GEMcon.git
cd GEMcon
```

### 2. Create necessary environments

```bash
conda env create -f envs/scanpy_env.yml
conda env create -f envs/cplex_env.yml
conda env create -f envs/gecko_env.yml
```

### 3. Usage
```bash
# Inside GEMcon directory
python scripts/cli.py --branch <transcriptomic|proteomic> [--task <annotate|metabolic>] [--algo <algorithm>] [--input <file>]
```
- `--branch` is required
- `--algo` is optional. Choices:
  - gimme 
  - tinit 
  - fastcore  (*testing*) 
  - imat      (*testing*) 
- `--input` is optional

### Optional

```bash
# Use as CLI instead
pip install -e . # install

gemcon --branch <transcriptomic|proteomic> [--task <annotate|metabolic>] [--algo <algorithm>] [--input <file>] # run!
```

### 4. Results

`results/` && `data_processed/`

---
## Status

In the works

### TODOs

GTex integration  
env. reproducibility  

## Docker

```bash
# pull the image
docker pull toluene123/gemcon:latest

docker run -it gemcon:latest bash # add --rm if you want post-exit deletion; -v for mounting a volume.

# check the config.yaml && ## run!
gemcon --help
```
