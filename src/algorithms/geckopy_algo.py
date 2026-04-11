import os
import argparse
import time
import requests
import pandas as pd
import geckopy
from geckopy.io import read_sbml_ec_model

from .base import BaseAlgorithm


class GeckopyAlgorithm(BaseAlgorithm):
    """Geckopy algorithm for enzyme-constrained metabolic model extraction"""

    def _parse_args(self):
        p = argparse.ArgumentParser()
        p.add_argument("--expr", help="Path to protein abundance file")
        p.add_argument("--model", help="Path to enzyme-constrained SBML model")
        p.add_argument("--format", help="Abundance file format: paxdb or generic", default=None)
        return p.parse_args()

    @property
    def algorithm_name(self):
        return "geckopy"

    # --- overrides: geckopy uses different config keys and file formats ---

    def _resolve_model_path(self):
        if self.args.model:
            return self.args.model.strip()
        return self.cfg.get("geckopy", {}).get("ec_model")

    def _detect_model_type(self):
        return "ec-model"

    def _resolve_expression_path(self):
        if self.args.expr:
            return self.args.expr.strip()
        return self.cfg.get("geckopy", {}).get("paxdb_data")

    def load_model(self):
        if not os.path.exists(self.model_path):
            raise FileNotFoundError(f"Model not found: {self.model_path}")
        return read_sbml_ec_model(self.model_path)

    def export_model(self, context_model, base_name):
        os.makedirs(self.output_dir, exist_ok=True)
        output_path = os.path.join(
            self.output_dir, f"geckopy_context_specific_model_{base_name}.xml"
        )
        geckopy.io.write_sbml_ec_model(context_model, output_path)
        return output_path

    # --- format detection ---

    def _detect_format(self):
        """Detect abundance file format from args, config, or file content"""
        fmt = getattr(self.args, "format", None) or self.cfg.get("geckopy", {}).get("format")
        if fmt:
            return fmt

        # Auto-detect: PaxDB files are tab-separated with # comments
        with open(self.expr_path) as f:
            first_lines = [f.readline() for _ in range(5)]
        if any(line.startswith("#") for line in first_lines):
            return "paxdb"
        return "generic"

    # --- format-specific loaders ---

    def _load_paxdb(self):
        """Load PaxDB tab-separated abundance file"""
        df = pd.read_csv(
            self.expr_path,
            sep="\t",
            comment="#",
            header=None,
            names=["internal_id", "string_external_id", "abundance"],
            usecols=[0, 1, 2],
        )
        # Strip STRING species prefix (e.g., 9606.ENSP00000370010 -> ENSP00000370010)
        df["ENSP"] = df["string_external_id"].str.replace("9606.", "", regex=False)
        return df

    def _map_to_uniprot(self, ensps):
        """Map Ensembl protein IDs to UniProt via ID mapping API"""
        try:
            job = requests.post(
                "https://rest.uniprot.org/idmapping/run",
                data={"from": "Ensembl_Protein", "to": "UniProtKB", "ids": ensps},
            ).json()
            job_id = job["jobId"]
        except Exception as e:
            raise RuntimeError(f"UniProt ID mapping failed: {e}") from e

        time.sleep(5)

        try:
            result = requests.get(
                f"https://rest.uniprot.org/idmapping/stream/{job_id}"
            ).json()
        except Exception as e:
            raise RuntimeError(f"UniProt results retrieval failed: {e}") from e

        map_df = pd.DataFrame(result.get("results", []))
        map_df = map_df.rename(columns={"from": "ENSP", "to": "UniProt"})
        return map_df

    def _prepare_protein_df(self, fmt):
        """Load abundance data and return DataFrame with protein_gecko_id, copies_per_cell, stdev"""
        if fmt == "paxdb":
            df = self._load_paxdb()
            ensps = df["ENSP"].drop_duplicates().tolist()

            # Map to UniProt
            map_df = self._map_to_uniprot(ensps)
            merged = df.merge(map_df, on="ENSP", how="left")
            print(f"Total rows mapped: {merged['UniProt'].notna().sum()}")

            prot = merged[["UniProt", "abundance"]].copy()
            prot = prot.dropna(subset=["UniProt"])
            prot = prot.rename(columns={"abundance": "copies_per_cell"})
            prot["protein_gecko_id"] = prot["UniProt"].apply(lambda x: f"prot_{x}[c]")
            prot["stdev"] = 0

        elif fmt == "generic":
            # Expects CSV with columns: uniprot, abundance (optional: stdev)
            df = pd.read_csv(self.expr_path)
            df.columns = df.columns.str.lower()

            if "uniprot" not in df.columns or "abundance" not in df.columns:
                raise ValueError("Generic format requires 'uniprot' and 'abundance' columns")

            prot = df[["uniprot", "abundance"]].copy()
            prot = prot.dropna(subset=["uniprot"])
            prot = prot.rename(columns={"uniprot": "UniProt", "abundance": "copies_per_cell"})
            prot["protein_gecko_id"] = prot["UniProt"].apply(lambda x: f"prot_{x}[c]")
            prot["stdev"] = df["stdev"] if "stdev" in df.columns else 0

        else:
            raise ValueError(f"Unknown format: {fmt}. Use 'paxdb' or 'generic'")

        return prot

    # --- run ---

    def run(self):
        # Load ec-model
        self.model = self.load_model()

        # Get geckopy config
        geckopy_cfg = self.cfg.get("geckopy", {})
        cell_params = geckopy_cfg.get("cell_params", {})
        vol = cell_params.get("vol", 2.3e-12)
        dens = cell_params.get("dens", 1.05)
        water = cell_params.get("water", 0.7)

        # Detect format and prepare protein data
        fmt = self._detect_format()
        print(f"Using abundance format: {fmt}")
        prot = self._prepare_protein_df(fmt)

        # Apply enzyme constraints
        ec_model_exp = geckopy.experimental.from_copy_number(
            self.model,
            index=prot["protein_gecko_id"],
            cell_copies=prot["copies_per_cell"],
            stdev=prot["stdev"],
            vol=vol,
            dens=dens,
            water=water,
        )

        print(f"Objective value: {ec_model_exp.slim_optimize()}")

        # Export
        base = os.path.splitext(os.path.basename(self.expr_path))[0]
        output_path = self.export_model(ec_model_exp, base)
        print(f"Saved to: {output_path}")

        return ec_model_exp