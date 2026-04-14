from abc import ABC, abstractmethod
import pandas as pd
import cobra
import yaml
import os


class BaseAlgorithm(ABC):
    """Base class for all model algorithms"""

    # from before
    MODEL_TYPES = {
        "human-gem": {"suffix": "_gencode.csv", "pattern": "__COBAMPGPRDOT__[0-9]{1}"},
        "recon3d": {"suffix": "_recon3d.csv", "pattern": "_AT[0-9]+$"},
    }

    def __init__(self, args=None, config_path="config.yaml"):
        self.args = args or self._parse_args()
        self.cfg = self._load_config(config_path)
        self.model_cfg = self.cfg.get("model_config", {})

        self.model_path = self._resolve_model_path()
        self.model_type = self._detect_model_type()
        self.expr_path = self._resolve_expression_path()
        self.output_dir = self._resolve_output_dir()

        self.model = None
        self.expression_data = None

    @abstractmethod
    def _parse_args(self):
        """Parse algorithm-specific command-line arguments"""
        pass

    @abstractmethod
    def run(self):
        """Execute the algorithm"""
        pass

    # common methods

    def _load_config(self, config_path):
        """Load YAML configuration file"""
        if os.path.exists(config_path):
            with open(config_path) as f:
                return yaml.safe_load(f) or {}
        return {}

    def _resolve_model_path(self):
        """Resolve model path from args or config"""
        if self.args.model:
            return self.args.model.strip()
        return self.cfg.get("model")

    def _detect_model_type(self):
        """Detect whether using Human-GEM or Recon3D generic model"""
        pattern = self.model_cfg.get(
            "alt_transcript_pattern", "__COBAMPGPRDOT__[0-9]{1}"
        )
        model_path = self.model_path or ""

        if (
            "_AT" in pattern
            or "recon3d" in model_path.lower()
            or "Recon3D" in model_path
        ):
            return "recon3d"
        return "human-gem"

    def _resolve_expression_path(self):
        """Resolve expression file path with auto-detection"""
        if self.args.expr:
            return self.args.expr.strip()

        suffix = self.MODEL_TYPES[self.model_type]["suffix"]
        data_dir = self.cfg.get("data_dir", "data/data_processed")

        files = [
            f
            for f in os.listdir(data_dir)
            if f.startswith("expression_data_") and f.endswith(suffix)
        ]

        if not files:
            return input(
                f"\nNo expression file found. Provide path ({self.model_type} format): "
            ).strip()

        paths = [os.path.join(data_dir, f) for f in files]
        return max(paths, key=os.path.getmtime)

    def _resolve_output_dir(self):
        """Resolve output directory from config"""
        return self.cfg.get("output", {}).get("models_dir", "./models")

    def load_model(self):
        """Load the metabolic model using Cobra for enzyme constrained ones"""
        if not os.path.exists(self.model_path):
            raise FileNotFoundError(f"Model not found: {self.model_path}")
        return cobra.io.read_sbml_model(self.model_path)

    def load_expression(self):
        """Load expression data as pandas df later switch to polars"""
        if not os.path.exists(self.expr_path):
            raise FileNotFoundError(f"Expression file not found: {self.expr_path}")
        return pd.read_csv(self.expr_path, index_col=0)

    def export_model(self, context_model, base_name):
        """Export context-specific model to SBML"""
        os.makedirs(self.output_dir, exist_ok=True)
        output_path = os.path.join(
            self.output_dir,
            f"{self.algorithm_name}_context_specific_model_{base_name}.xml",
        )
        cobra.io.write_sbml_model(context_model, output_path)
        return output_path

    def create_context_model(self, model, selected_reaction_ids):
        """Create context-specific model by removing unselected reactions as per troppo"""
        ctx_model = model.copy()
        selected_ids = set(selected_reaction_ids)
        to_remove = [r for r in ctx_model.reactions if r.id not in selected_ids]
        ctx_model.remove_reactions(to_remove, remove_orphans=True)
        return ctx_model

    @property
    @abstractmethod
    def algorithm_name(self):
        """Return algorithm name for file naming"""
        pass
