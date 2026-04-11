import pandas as pd
import numpy as np
import cobra
import re
import os
from cobamp.core.linear_systems import get_default_solver
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper, ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.imat import IMAT, IMATProperties
import optlang
from optlang.cplex_interface import Model as CPLEXModel

print("COBAMP default solver:", get_default_solver())

# parsing rule 
patt = re.compile('__COBAMPGPRDOT__[0-9]{1}') # e.g __COBAMPGPRDOT__2
replace_alt_transcripts = lambda x: patt.sub('', x) #  empty string (prune)

# load model and expression data
model = cobra.io.read_sbml_model('/home/biodata/aman/Human-GEM/model/Human-GEM.xml')

expression_data = pd.read_csv('/home/biodata/aman/data/data_processed/Brain_fibroblasts_gencode.csv', index_col=0)

expression_data_transposed = expression_data.T

print("Transposed:", expression_data_transposed.head())

omics_container = TabularReader(path_or_df=expression_data_transposed, 
                                nomenclature='gene',
                                omics_type='transcriptomics').to_containers()

single_sample = omics_container[0]

model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999,
                                      gpr_gene_parse_function=replace_alt_transcripts)

data_map = single_sample.get_integrated_data_map(model_reader=model_wrapper.model_reader,
                                                 and_func=min, or_func=sum)

def score_apply(reaction_map_scores):
    return {k:0  if v is None else v for k, v in reaction_map_scores.items()}

continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)

scores = continuous_integration.integrate(data_map=data_map)

reaction_ids = model_wrapper.model_reader.r_ids

exp_vector = np.array([scores[rid] for rid in reaction_ids])

# redundant
os.environ["COBAMP_SOLVER"] = "CPLEX"
optlang.config = {'solver': 'cplex'}
os.environ["OPTLANG_DEFAULT_SOLVER"] = "CPLEX"

# Create the properties for the IMAT algorithm.
properties = IMATProperties(exp_vector=exp_vector, exp_thresholds=(25,75))

# Run the iMAT algorithm.
imat = IMAT(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)

print('Running iMAT')
model_imat = imat.run()
print("\nReaction indices to keep: ", pd.DataFrame(model_imat))

pd.DataFrame(model_imat).to_csv('troppo_output_imat.tsv', sep='\t')

selected_reactions = [model.reactions[i] for i in model_imat.flatten().tolist()]

print("\nSelected reactions: ", selected_reactions)

ctx_model = model.copy()

selected_ids = [r.id for r in selected_reactions]

to_remove = [r for r in ctx_model.reactions if r.id not in selected_ids]

ctx_model.remove_reactions(to_remove, remove_orphans=True)

print(ctx_model)

print(len(model.reactions), len(ctx_model.reactions))
# export
cobra.io.write_sbml_model(ctx_model, "imat_context_specific_model.xml")