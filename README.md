# GEMcon

[![test_ga](https://github.com/ToruOkadaOi/GEMcon/actions/workflows/test.yml/badge.svg)](https://github.com/ToruOkadaOi/GEMcon/actions/workflows/test.yml) [![Docker](https://img.shields.io/docker/v/toluene123/gemcon?label=toluene123/gemcon)](https://hub.docker.com/r/toluene123/gemcon)

Build **context-specific human metabolic models** from expression and proteomic abundance data, based on [Human-GEM](https://github.com/SysBioChalmers/Human-GEM).

GEMcon does data fetching, processing, and model building using multiple reconstruction algorithms (GIMME, tINIT, RIPTiDe, CORDA, and more) under one roof!

---

## Data Sources

GEMcon integrates abundance data from:

- **scRNA-seq** — [Human Cell Atlas](https://www.humancellatlas.org)
- **Bulk RNA-seq** — [GTEx](https://gtexportal.org/home/) # check utils
- **Proteomics** — [PaxDB](https://pax-db.org)

---

## Prerequisites

- Python 3.x
- [Conda](https://docs.conda.io/en/latest/) (Miniconda or Anaconda)
- A valid [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) license (required by several reconstruction algorithms)

---

## Installation

### Option A: Docker (recommended)

The Docker image ships with all dependencies pre-configured — no conda or CPLEX setup required.
```bash
docker pull toluene123/gemcon:latest

docker run -it toluene123/gemcon:latest bash
# add --rm for post-exit cleanup; -v to mount a volume for results

gemcon --help
```

### Option B: Install from source
```bash
git clone https://github.com/ToruOkadaOi/GEMcon.git
cd GEMcon

# Create the required conda environments
conda env create -f envs/scanpy_env.yml
conda env create -f envs/cplex_env.yml
conda env create -f envs/gecko_env.yml

# (Optional) Install as a CLI tool
pip install -e .
```

---

## Usage
```bash
# If installed via pip:
gemcon --branch <transcriptomic|proteomic> [--task <annotate|metabolic>] [--algo <algorithm>] [--input <file>]

# Or run directly from the repo:
python scripts/cli.py --branch <transcriptomic|proteomic> [--task <annotate|metabolic>] [--algo <algorithm>] [--input <file>]
```

### Arguments

| Argument    | Required | Description |
|-------------|----------|-------------|
| `--branch`  | Yes      | Data type to use: `transcriptomic` or `proteomic` |
| `--task`    | No       | Pipeline stage: `annotate` (process expression data) or `metabolic` (build model) |
| `--algo`    | No       | Reconstruction algorithm (see below) |
| `--input`   | No       | Path to a custom input file |

### Supported Algorithms

| Algorithm  | Status  |
|------------|---------|
| `gimme`    | Stable  |
| `tinit`    | Stable  |
| `riptide`  | Stable  |
| `corda`    | Stable  |
| `fastcore` | Testing |
| `imat`     | Testing |

---

## Output

Results are written to two directories:

- **`results/`** — Final context-specific metabolic models
- **`data_processed/`** — Intermediate processed expression/abundance data

---

## Roadmap


See the [Issues](https://github.com/ToruOkadaOi/GEMcon/issues) page for more details.

---

> **Status:** This project is no longer actively maintained.

## License

This project is licensed under the [GPL-3.0 License](LICENSE).
