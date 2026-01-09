import pandas as pd
import numpy as np
import cobra
import re
from pathlib import Path

# Troppo imports
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.gimme import GIMME as TroppoGIMME, GIMMEProperties  # Keep this line

from .base import AlgorithmBase

class GIMME(AlgorithmBase):    
    def __init__(self, config):
        super().__init__(config)
        pattern = self.config.get('model_config', {}).get(
        'alt_transcript_pattern', 
        '__COBAMPGPRDOT__[0-9]{1}'
        )
        self.patt = re.compile(pattern)
        self.replace_alt_transcripts = lambda x: self.patt.sub('', x)

    
    def run(self, model, expression_data):
        expression_data_T = expression_data.T

        omics_container = TabularReader(
        path_or_df=expression_data_T,
        nomenclature='gene',
        omics_type='transcriptomics'
        ).to_containers()

        single_sample = omics_container[0]

        model_wrapper = ReconstructionWrapper(
            model=model,
            ttg_ratio=9999,
            gpr_gene_parse_function=self.replace_alt_transcripts
        )

        data_map = single_sample.get_integrated_data_map(
            model_reader=model_wrapper.model_reader,
            and_func=min,
            or_func=sum
        )

        def score_apply(reaction_map_scores):
            return {k: 0 if v is None else v for k, v in reaction_map_scores.items()}
        
        continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
        scores = continuous_integration.integrate(data_map=data_map)
        
        gimme_cfg = self.config.get('gimme', {})
        obj_frac = gimme_cfg.get('obj_frac', 0.1)
        flux_threshold = gimme_cfg.get('flux_threshold', 0.25)
        solver = gimme_cfg.get('solver', 'CPLEX')
        
        objective_rxn = self.config.get('model_config', {}).get('objective_reaction')
        if not objective_rxn:
            raise ValueError("objective_reaction must be specified in config")      

        idx_objective = model_wrapper.model_reader.r_ids.index(objective_rxn)
        
        properties = GIMMEProperties(
            exp_vector=[v for k, v in scores.items()],
            obj_frac=obj_frac,
            objectives=[{idx_objective: 1}],
            preprocess=True,
            flux_threshold=flux_threshold,
            solver=solver,
            reaction_ids=model_wrapper.model_reader.r_ids,
            metabolite_ids=model_wrapper.model_reader.m_ids
        )
        
        gimme_algo = TroppoGIMME(
            S=model_wrapper.S,
            lb=model_wrapper.lb,
            ub=model_wrapper.ub,
            properties=properties
        )
        gimme_result = gimme_algo.run()

        selected_reactions = [model.reactions[i] for i in gimme_result]
        ctx_model = model.copy()

        selected_ids = [r.id for r in selected_reactions]
        to_remove = [r for r in ctx_model.reactions if r.id not in selected_ids]
        
        ctx_model.remove_reactions(to_remove, remove_orphans=True)
        
        self.result = gimme_result
        self.model = ctx_model
        
        return ctx_model