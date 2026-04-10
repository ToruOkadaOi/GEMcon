import os
import csv
import argparse
import pandas as pd
from corda import CORDA, reaction_confidence

from .base import BaseAlgorithm


class CordaAlgorithm(BaseAlgorithm):
    """CORDA algorithm for context-specific metabolic model extraction"""
    
    def _parse_args(self):
        p = argparse.ArgumentParser()
        p.add_argument("--expr", help="Path to expression data csv")
        p.add_argument("--model", help="Path to SBML model")
        return p.parse_args()
    
    @property
    def algorithm_name(self):
        return "corda"
    
    def run(self):
        # Load model and expression
        self.model = self.load_model()
        self.model.solver = "cplex"  # CORDA also works with glpk
        self.expression_data = self.load_expression()
        
        # Get CORDA config
        corda_cfg = self.cfg.get("corda", {})
        
        # Get expression values (single sample or mean across samples)
        if self.expression_data.shape[1] == 1:
            gene_expr = self.expression_data.iloc[:, 0]
        else:
            gene_expr = self.expression_data.mean(axis=1)
        
        # Confidence thresholds — use absolute if provided, otherwise compute from quantiles
        threshold_absent = corda_cfg.get(
            "threshold_absent", gene_expr.quantile(corda_cfg.get("quantile_absent", 0.05))
        )
        threshold_unknown = corda_cfg.get(
            "threshold_unknown", gene_expr.quantile(corda_cfg.get("quantile_unknown", 0.25))
        )
        threshold_low = corda_cfg.get(
            "threshold_low", gene_expr.quantile(corda_cfg.get("quantile_low", 0.50))
        )
        threshold_medium = corda_cfg.get(
            "threshold_medium", gene_expr.quantile(corda_cfg.get("quantile_medium", 0.75))
        )
        
        # Map gene expression to confidence levels (-1 to 3)
        def expr_to_confidence(val):
            if val <= threshold_absent:
                return -1  # absent
            elif val <= threshold_unknown:
                return 0   # unknown
            elif val <= threshold_low:
                return 1   # low confidence
            elif val <= threshold_medium:
                return 2   # medium confidence
            else:
                return 3   # high confidence
        
        # Build gene confidence dictionary
        conf_genes = {}
        for gene_id, expr_val in gene_expr.items():
            # Strip transcript version if present (e.g., ENSG000001234.5 -> ENSG000001234)
            gene_clean = str(gene_id).split(".")[0]
            conf_genes[gene_clean] = expr_to_confidence(expr_val)
        
        # Map gene confidences to reaction confidences using GPR rules
        conf_reactions = {}
        for rxn in self.model.reactions:
            if rxn.gene_reaction_rule:
                conf_reactions[rxn.id] = reaction_confidence(rxn, conf_genes)
            else:
                # No GPR rule — unknown confidence
                conf_reactions[rxn.id] = 0
        
        # Run CORDA
        met_prod = corda_cfg.get("met_prod", [])
        corda_obj = CORDA(self.model, conf_reactions, met_prod=met_prod)
        corda_obj.build()
        
        # Build context-specific model
        ctx_rxns = [r.id for r in self.model.reactions if corda_obj.included.get(r.id, False)]
        ctx_model = self.create_context_model(self.model, ctx_rxns)
        
        print(f"Base model reactions: {len(self.model.reactions)}")
        print(f"Context model reactions: {len(ctx_model.reactions)}")
        
        # Export model
        base = os.path.splitext(os.path.basename(self.expr_path))[0]
        output_path = self.export_model(ctx_model, base)
        print(f"Saved to: {output_path}")
        
        # Save reaction confidences
        conf_path = os.path.join(self.output_dir, f"corda_reaction_confidences_{base}.tsv")
        with open(conf_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["reaction_id", "confidence"])
            for k, v in conf_reactions.items():
                writer.writerow([k, v])
        print(f"Saved confidences to: {conf_path}")
        
        solution = ctx_model.optimize()
        print(f"Objective value: {solution.objective_value}")
        
        return ctx_model