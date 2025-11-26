# GEMcon

Pipeline for fetching expression~~/proteomic abundance~~ data and building **context-specific human metabolic models** (Human-GEM  based)

---

## Overview

GEMcon integrates data from:
- scRNA-seq (Human Cell Atlas)
- bulk RNA-seq (GTEx) -->proposed
- proteomics (PaxDB) -->proposed

---

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/aman/GEMcon
cd GEMcon
```

### 2. Create necessary environments

```bash
conda env create -f envs/scanpy_legacy.yml
conda env create -f envs/cplex.yml
conda env create -f envs/gecko.yml
```

### 3. Usage
```bash
# Inside GEMcon directory
python scripts/run.py --branch <celltype_annotated|metabolic> [--input <file>]
```
- `--branch` is required
- `--input` is optional

### 4. Results

`results/` && `data_processed/`

---
## Status

Early stage

### TODOs
