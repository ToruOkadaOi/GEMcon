"""Build a toy COBRA model for testing algorithms."""

import cobra
import pandas as pd
import os


def build_toy_model():
    """Create a minimal COBRA model with 5 reactions, 4 metabolites, 3 genes."""
    model = cobra.Model("toy_model")

    m1 = cobra.Metabolite("M1_c", compartment="c")
    m2 = cobra.Metabolite("M2_c", compartment="c")
    m3 = cobra.Metabolite("M3_c", compartment="c")
    m4 = cobra.Metabolite("M4_c", compartment="c")

    r1 = cobra.Reaction("R1", name="M1 -> M2")
    r1.add_metabolites({m1: -1, m2: 1})
    r1.lower_bound = 0
    r1.upper_bound = 1000

    r2 = cobra.Reaction("R2", name="M2 -> M3")
    r2.add_metabolites({m2: -1, m3: 1})
    r2.lower_bound = 0
    r2.upper_bound = 1000

    r3 = cobra.Reaction("R3", name="M3 -> M4")
    r3.add_metabolites({m3: -1, m4: 1})
    r3.lower_bound = 0
    r3.upper_bound = 1000

    r4 = cobra.Reaction("R4", name="M4 -> M1")
    r4.add_metabolites({m4: -1, m1: 1})
    r4.lower_bound = -1000
    r4.upper_bound = 1000

    r5 = cobra.Reaction("R5", name="sink M3")
    r5.add_metabolites({m3: -1})
    r5.lower_bound = 0
    r5.upper_bound = 1000

    model.add_reactions([r1, r2, r3, r4, r5])
    model.add_metabolites([m1, m2, m3, m4])

    r1.gene_reaction_rule = "G1"
    r2.gene_reaction_rule = "G1 and G2"
    r3.gene_reaction_rule = "G2"
    r4.gene_reaction_rule = "G3"
    r5.gene_reaction_rule = "G1 or G3"

    model.objective = "R1"

    return model


def build_toy_expression():
    """Create matching expression CSV with 3 genes."""
    expr_data = pd.DataFrame(
        {"sample_1": [10.0, 5.0, 2.0]}, index=pd.Index(["G1", "G2", "G3"], name="gene")
    )
    return expr_data


def main():
    model = build_toy_model()
    model_path = os.path.join(os.path.dirname(__file__), "toy_model.xml")
    cobra.io.write_sbml_model(model, model_path)
    print(f"Saved model to {model_path} with {len(model.reactions)} reactions")

    expr = build_toy_expression()
    expr_path = os.path.join(os.path.dirname(__file__), "toy_expr.csv")
    expr.to_csv(expr_path)
    print(f"Saved expression to {expr_path}")


if __name__ == "__main__":
    main()
