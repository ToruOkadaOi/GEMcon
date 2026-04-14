import re
import os
import argparse
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties

from .base import BaseAlgorithm


class TinitAlgorithm(BaseAlgorithm):
    """tINIT algorithm for context-specific metabolic model extraction"""

    def _parse_args(self):
        p = argparse.ArgumentParser()
        p.add_argument("--expr", help="Path to expression data csv")
        p.add_argument("--model", help="Path to SBML model")
        return p.parse_args()

    @property
    def algorithm_name(self):
        return "tinit"

    def run(self):
        self.model = self.load_model()
        self.expression_data = self.load_expression()

        pattern = self.model_cfg.get(
            "alt_transcript_pattern", "__COBAMPGPRDOT__[0-9]{1}"
        )
        patt = re.compile(pattern)
        replace_alt_transcripts = lambda x: patt.sub("", x)

        expr_transposed = self.expression_data.T

        tinit_cfg = self.cfg.get("tinit", {})
        nomenclature = tinit_cfg.get("nomenclature", "ensembl_gene_id")

        omics_container = TabularReader(
            path_or_df=expr_transposed,
            nomenclature=nomenclature,
            omics_type="transcriptomics",
        ).to_containers()

        single_sample = omics_container[0]

        model_wrapper = ReconstructionWrapper(
            model=self.model,
            ttg_ratio=9999,
            gpr_gene_parse_function=replace_alt_transcripts,
        )

        data_map = single_sample.get_integrated_data_map(
            model_reader=model_wrapper.model_reader, and_func=min, or_func=sum
        )

        def score_apply(reaction_map_scores):
            return {k: 0 if v is None else v for k, v in reaction_map_scores.items()}

        continuous_integration = ContinuousScoreIntegrationStrategy(
            score_apply=score_apply
        )
        scores = continuous_integration.integrate(data_map=data_map)

        # Get tINIT-specific parameters
        tinit_cfg = self.cfg.get("tinit", {})

        essential_rxn_ids = tinit_cfg.get("essential_reactions", [])
        if not essential_rxn_ids:
            user_input = input("\nEnter essential reaction ids: ").strip()
            essential_rxn_ids = [r.strip() for r in user_input.split(",") if r.strip()]

        essential_reactions = [
            model_wrapper.model_reader.r_ids.index(rid) for rid in essential_rxn_ids
        ]

        solver = tinit_cfg.get("solver", "CPLEX")

        properties = tINITProperties(
            reactions_scores=[v for k, v in scores.items()],
            solver=solver,
            essential_reactions=essential_reactions,
        )

        tinit = tINIT(
            S=model_wrapper.S,
            lb=model_wrapper.lb,
            ub=model_wrapper.ub,
            properties=properties,
        )

        model_tinit = tinit.run()

        base = os.path.splitext(os.path.basename(self.expr_path))[0]

        selected_reactions = [
            self.model.reactions[i] for i in model_tinit.flatten().tolist()
        ]
        selected_ids = [r.id for r in selected_reactions]

        ctx_model = self.create_context_model(self.model, selected_ids)

        output_path = self.export_model(ctx_model, base)
        print(f"Saved to: {output_path}")

        solution = ctx_model.optimize()
        print(f"Objective value: {solution.objective_value}")

        return ctx_model
