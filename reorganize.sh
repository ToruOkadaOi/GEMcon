#!/bin/bash

mkdir -p data_raw data_processed models scripts notebooks envs external_tools results/figures results/logs

# RAW DATA
mv expression_data.csv data_raw/
mv gencode.v43.annotation.gtf data_raw/
mv gencode.v49.annotation.gtf.gz data_raw/
mv t2.loom data_raw/
mv HCA_downloads data_raw/
mv misc data_raw/
mv data data_raw/
mv master.zip data_raw/

# PROCESSED DATA
mv expression_by_celltype data_processed/
mv expression_data_ensembl.csv data_processed/
mv gtf_reference_genesymbols_ensemblid.tsv data_processed/
mv processed data_processed/
mv organized data_processed/

# MODELS
mv fastcore_context_specific_model.xml models/
mv fastcore_output_tinit.tsv models/
mv gimme_context_specific_model.xml models/
mv gimme_output_tinit.tsv models/
mv tinit_context_specific_model.xml models/
mv troppo_output_tinit.tsv models/
mv tINIT_test.lp models/
# clean duplicates from scripts
mv scripts/gimme_context_specific_model.xml models/ 2>/dev/null
mv scripts/gimme_output_tinit.tsv models/ 2>/dev/null
mv scripts/tinit_context_specific_model.xml models/ 2>/dev/null
mv scripts/tINIT_test.lp models/ 2>/dev/null

# NOTEBOOKS
mv scripts/*.ipynb notebooks/

# SCRIPTS
mv scripts/*.py scripts/
mv scripts/*.sh scripts/

# ENV FILES
mv scanpy_legacy.yml envs/

# EXTERNAL TOOLS
mv cobratoolbox-master external_tools/
