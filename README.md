# GEMcon

[![test_ga](https://github.com/ToruOkadaOi/GEMcon/actions/workflows/test.yml/badge.svg)](https://github.com/ToruOkadaOi/GEMcon/actions/workflows/test.yml) [![Docker](https://img.shields.io/docker/v/toluene123/gemcon?label=toluene123/gemcon)](https://hub.docker.com/r/toluene123/gemcon)

Build context-specific human metabolic models from expression and protein abundance data. Wraps several reconstruction algorithms (GIMME, tINIT, RIPTiDe, CORDA, FASTCORE, iMAT) on top of [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) with a common CLI for data fetching, normalization, and model extraction.

## Data sources

- scRNA-seq — [Human Cell Atlas](https://www.humancellatlas.org)
- Bulk RNA-seq — [GTEx](https://gtexportal.org/home/)
- Proteomics — [PaxDB](https://pax-db.org)

## Requirements

- Python 3.10
- Conda (Miniconda or Anaconda)
- CPLEX 22.1+ with a valid license (required by some algorithms)

## Installation

### Docker

The image bundles the runtime. CPLEX is not included; mount your license and run `install_cplex.sh` after starting the container.

```bash
docker pull toluene123/gemcon:latest
docker run -it -v $PWD/results:/app/results toluene123/gemcon:latest bash
gemcon --help
```

### From source

```bash
git clone https://github.com/ToruOkadaOi/GEMcon.git
cd GEMcon

conda env create -f envs/environment.yml
conda activate gemcon

# if pip complains about numpy:
pip install --no-deps --force-reinstall troppo==0.0.7 cobamp==0.2.1

# CPLEX env vars (adjust the path to your install). Add to ~/.bashrc to persist.
export CPLEX_STUDIO_DIR2212=/path/to/CPLEX_Studio2212
export LD_LIBRARY_PATH=$CPLEX_STUDIO_DIR2212/cplex/bin/x86-64_linux:$LD_LIBRARY_PATH

# Optional: install as a CLI
pip install -e .
```

## Usage

```bash
gemcon --branch <transcriptomic|proteomic> [--task <annotate|metabolic>] [--algo <algorithm>] [--input <file>]

# Or directly:
python -m scripts.cli --branch transcriptomic --task metabolic --algo gimme
```

| Argument   | Required | Description |
|------------|----------|-------------|
| `--branch` | yes      | `transcriptomic` or `proteomic` |
| `--task`   | no       | `annotate` or `metabolic` |
| `--algo`   | no       | Reconstruction algorithm |
| `--input`  | no       | Path to input file (otherwise resolved from config) |

### Algorithms

| Algorithm  | Status |
|------------|--------|
| `gimme`    | working |
| `fastcore` | working |
| `tinit`    | working |
| `riptide`  | working |
| `corda`    | working |
| `geckopy`  | working but requires a populated ec-model |
| `imat`     | testing |

## Output

- `models/` — context-specific SBML files
- `data/data_processed/` — intermediate normalized/pooled expression data

## License

[GPL-3.0](LICENSE)