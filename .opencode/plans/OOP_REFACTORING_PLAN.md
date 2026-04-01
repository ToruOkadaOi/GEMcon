# GEMcon OOP Refactoring Plan

## Overview

This document outlines a step-by-step plan to refactor the algorithm scripts in `src/` from procedural code to an object-oriented structure using inheritance. The goal is to eliminate code duplication (~70% of each file is identical) while maintaining flexibility for algorithm-specific behavior.

**Estimated Total Time: 6-8 hours**

---

## Current State Analysis

### Files to Refactor
| File | Lines | Algorithm |
|------|-------|------------|
| `src/gimme.py` | 184 | GIMME |
| `src/tinit.py` | 147 | tINIT |
| `src/corda_algo.py` | 232 | CORDA |
| `src/fastcore_beta.py` | 82 | FastCORE |
| `src/riptide_algo.py` | ~200+ | RIPTiDe |

### Code Duplication Pattern

Each algorithm file contains nearly identical:

```python
# 1. CONFIG LOADING (~15 lines)
cfg = {}
if os.path.exists("config.yaml"):
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f) or {}

# 2. MODEL PATH RESOLUTION (~20 lines)
if args.model:
    model_path = args.model.strip()
else:
    model_path = cfg.get("model")

# Detect model type (Human-GEM vs Recon3D)
pattern = model_cfg.get('alt_transcript_pattern')
# ... pattern matching logic ...

# 3. EXPRESSION FILE LOOKUP (~20 lines)
if args.expr:
    expr_path = args.expr.strip()
else:
    files = [f for f in os.listdir("data/data_processed") 
             if f.startswith("expression_data_") and f.endswith(expr_suffix)]
    # ... auto-detection logic ...

# 4. MODEL LOADING (~5 lines)
model = cobra.io.read_sbml_model(model_path)

# 5. EXPRESSION LOADING (~5 lines)
expression_data = pd.read_csv(expr_path, index_col=0)

# 6. OUTPUT DIRECTORY CREATION (~5 lines)
mod_dir = cfg.get('output', {}).get('models_dir', './models')
os.makedirs(mod_dir, exist_ok=True)

# 7. MODEL EXPORT (~3 lines)
cobra.io.write_sbml_model(ctx_model, output_path)
```

**Total duplicate code: ~120 lines per file**

---

## Target Architecture

```
src/algorithms/
├── __init__.py           # Exports all algorithm classes
├── base.py               # BaseAlgorithm abstract class
├── gimme.py              # GimmeAlgorithm(BaseAlgorithm)
├── tinit.py              # TinitAlgorithm(BaseAlgorithm)
├── corda.py              # CordaAlgorithm(BaseAlgorithm)
├── fastcore.py           # FastcoreAlgorithm(BaseAlgorithm)
└── riptide.py            # RiptideAlgorithm(BaseAlgorithm)
```

### BaseAlgorithm Class

```python
# src/algorithms/base.py
from abc import ABC, abstractmethod
import pandas as pd
import cobra
import yaml
import os
import argparse
from pathlib import Path

class BaseAlgorithm(ABC):
    """Base class for all context-specific metabolic model algorithms."""
    
    MODEL_TYPES = {
        'human-gem': {'suffix': '_gencode.csv', 'pattern': '__COBAMPGPRDOT__[0-9]{1}'},
        'recon3d': {'suffix': '_recon3d.csv', 'pattern': '_AT[0-9]+$'}
    }
    
    def __init__(self, args=None, config_path="config.yaml"):
        self.args = args or self._parse_args()
        self.cfg = self._load_config(config_path)
        self.model_cfg = self.cfg.get('model_config', {})
        
        self.model_path = self._resolve_model_path()
        self.model_type = self._detect_model_type()
        self.expr_path = self._resolve_expression_path()
        self.output_dir = self._resolve_output_dir()
        
        self.model = None
        self.expression_data = None
    
    @abstractmethod
    def _parse_args(self):
        """Parse algorithm-specific command-line arguments."""
        pass
    
    @abstractmethod
    def run(self):
        """Execute the algorithm. Must be implemented by subclasses."""
        pass
    
    # === COMMON METHODS ===
    
    def _load_config(self, config_path):
        """Load YAML configuration file."""
        if os.path.exists(config_path):
            with open(config_path) as f:
                return yaml.safe_load(f) or {}
        return {}
    
    def _resolve_model_path(self):
        """Resolve model path from args or config."""
        if self.args.model:
            return self.args.model.strip()
        return self.cfg.get('model')
    
    def _detect_model_type(self):
        """Detect whether using Human-GEM or Recon3D."""
        pattern = self.model_cfg.get('alt_transcript_pattern', '__COBAMPGPRDOT__[0-9]{1}')
        model_path = self.model_path or ''
        
        if '_AT' in pattern or 'recon3d' in model_path.lower() or 'Recon3D' in model_path:
            return 'recon3d'
        return 'human-gem'
    
    def _resolve_expression_path(self):
        """Resolve expression file path with auto-detection."""
        if self.args.expr:
            return self.args.expr.strip()
        
        suffix = self.MODEL_TYPES[self.model_type]['suffix']
        data_dir = self.cfg.get('data_dir', 'data/data_processed')
        
        files = [f for f in os.listdir(data_dir) 
                 if f.startswith("expression_data_") and f.endswith(suffix)]
        
        if not files:
            return input(f"\nNo expression file found. Provide path ({self.model_type} format): ").strip()
        
        paths = [os.path.join(data_dir, f) for f in files]
        return max(paths, key=os.path.getmtime)
    
    def _resolve_output_dir(self):
        """Resolve output directory from config."""
        return self.cfg.get('output', {}).get('models_dir', './models')
    
    def load_model(self):
        """Load the metabolic model using Cobra."""
        if not os.path.exists(self.model_path):
            raise FileNotFoundError(f"Model not found: {self.model_path}")
        return cobra.io.read_sbml_model(self.model_path)
    
    def load_expression(self):
        """Load expression data as DataFrame."""
        if not os.path.exists(self.expr_path):
            raise FileNotFoundError(f"Expression file not found: {self.expr_path}")
        return pd.read_csv(self.expr_path, index_col=0)
    
    def export_model(self, context_model, base_name):
        """Export context-specific model to SBML."""
        os.makedirs(self.output_dir, exist_ok=True)
        output_path = os.path.join(self.output_dir, f"{self.algorithm_name}_context_specific_model_{base_name}.xml")
        cobra.io.write_sbml_model(context_model, output_path)
        return output_path
    
    def create_context_model(self, model, selected_reaction_ids):
        """Create context-specific model by removing unselected reactions."""
        ctx_model = model.copy()
        selected_ids = set(selected_reaction_ids)
        to_remove = [r for r in ctx_model.reactions if r.id not in selected_ids]
        ctx_model.remove_reactions(to_remove, remove_orphans=True)
        return ctx_model
    
    @property
    @abstractmethod
    def algorithm_name(self):
        """Return algorithm name for file naming."""
        pass
```

### GimmeAlgorithm Implementation

```python
# src/algorithms/gimme.py
import re
import os
import argparse
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties

from .base import BaseAlgorithm

class GimmeAlgorithm(BaseAlgorithm):
    """GIMME algorithm for context-specific metabolic model extraction."""
    
    def _parse_args(self):
        p = argparse.ArgumentParser()
        p.add_argument("--expr", help="Path to expression data csv")
        p.add_argument("--model", help="Path to SBML model")
        return p.parse_args()
    
    @property
    def algorithm_name(self):
        return "gimme"
    
    def run(self):
        # Load model and expression
        self.model = self.load_model()
        self.expression_data = self.load_expression()
        
        # Get gene pattern for GPR parsing
        pattern = self.model_cfg.get('alt_transcript_pattern', '__COBAMPGPRDOT__[0-9]{1}')
        patt = re.compile(pattern)
        replace_alt_transcripts = lambda x: patt.sub('', x)
        
        # Prepare expression data for Troppo
        expr_transposed = self.expression_data.T
        
        omics_container = TabularReader(
            path_or_df=expr_transposed,
            nomenclature='gene',
            omics_type='transcriptomics'
        ).to_containers()
        
        single_sample = omics_container[0]
        
        model_wrapper = ReconstructionWrapper(
            model=self.model,
            ttg_ratio=9999,
            gpr_gene_parse_function=replace_alt_transcripts
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
        
        # Get algorithm-specific parameters
        gimme_cfg = self.cfg.get('gimme', {})
        
        obj_frac = gimme_cfg.get('obj_frac', 0.1)
        flux_threshold = gimme_cfg.get('flux_threshold', 0.25)
        solver = gimme_cfg.get('solver', 'CPLEX')
        
        # Get objective reaction
        objective_rxn = self.model_cfg.get('objective_reaction', 'MAR02388')
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
        
        # Run GIMME
        gimme = GIMME(
            S=model_wrapper.S,
            lb=model_wrapper.lb,
            ub=model_wrapper.ub,
            properties=properties
        )
        
        gimme_run = gimme.run()
        
        # Create and export context model
        base = os.path.splitext(os.path.basename(self.expr_path))[0]
        
        selected_reactions = [self.model.reactions[i] for i in gimme_run]
        selected_ids = [r.id for r in selected_reactions]
        
        ctx_model = self.create_context_model(self.model, selected_ids)
        
        print(f"Base model reactions: {len(self.model.reactions)}")
        print(f"Context model reactions: {len(ctx_model.reactions)}")
        
        output_path = self.export_model(ctx_model, base)
        print(f"Saved to: {output_path}")
        
        solution = ctx_model.optimize()
        print(f"Objective value: {solution.objective_value}")
        
        return ctx_model
```

### tINIT Implementation

```python
# src/algorithms/tinit.py
import re
import os
import argparse
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties

from .base import BaseAlgorithm

class TinitAlgorithm(BaseAlgorithm):
    """tINIT algorithm for context-specific metabolic model extraction."""
    
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
        
        pattern = self.model_cfg.get('alt_transcript_pattern', '__COBAMPGPRDOT__[0-9]{1}')
        patt = re.compile(pattern)
        replace_alt_transcripts = lambda x: patt.sub('', x)
        
        expr_transposed = self.expression_data.T
        
        tinit_cfg = self.cfg.get('tinit', {})
        nomenclature = tinit_cfg.get('nomenclature', 'ensembl_gene_id')
        
        omics_container = TabularReader(
            path_or_df=expr_transposed,
            nomenclature=nomenclature,
            omics_type='transcriptomics'
        ).to_containers()
        
        single_sample = omics_container[0]
        
        model_wrapper = ReconstructionWrapper(
            model=self.model,
            ttg_ratio=9999,
            gpr_gene_parse_function=replace_alt_transcripts
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
        
        # Get tINIT-specific parameters
        tinit_cfg = self.cfg.get('tinit', {})
        
        essential_rxn_ids = tinit_cfg.get('essential_reactions', [])
        if not essential_rxn_ids:
            user_input = input("\nEnter essential reaction ids: ").strip()
            essential_rxn_ids = [r.strip() for r in user_input.split(',') if r.strip()]
        
        essential_reactions = [model_wrapper.model_reader.r_ids.index(rid) for rid in essential_rxn_ids]
        
        solver = tinit_cfg.get('solver', 'CPLEX')
        
        properties = tINITProperties(
            reactions_scores=[v for k, v in scores.items()],
            solver=solver,
            essential_reactions=essential_reactions
        )
        
        tinit = tINIT(
            S=model_wrapper.S,
            lb=model_wrapper.lb,
            ub=model_wrapper.ub,
            properties=properties
        )
        
        model_tinit = tinit.run()
        
        base = os.path.splitext(os.path.basename(self.expr_path))[0]
        
        selected_reactions = [self.model.reactions[i] for i in model_tinit.flatten().tolist()]
        selected_ids = [r.id for r in selected_reactions]
        
        ctx_model = self.create_context_model(self.model, selected_ids)
        
        output_path = self.export_model(ctx_model, base)
        print(f"Saved to: {output_path}")
        
        solution = ctx_model.optimize()
        print(f"Objective value: {solution.objective_value}")
        
        return ctx_model
```

---

## Step-by-Step Implementation Plan

### Phase 1: Create Base Class (1 hour)
1. Create `src/algorithms/` directory
2. Create `src/algorithms/__init__.py`
3. Create `src/algorithms/base.py` with `BaseAlgorithm` class
4. Test base class loads without errors

### Phase 2: Refactor GIMME (1.5 hours)
1. Create `src/algorithms/gimme.py`
2. Implement `GimmeAlgorithm(BaseAlgorithm)`
3. Move algorithm-specific logic from `src/gimme.py`
4. Test: run both old and new, compare outputs

### Phase 3: Refactor tINIT (1.5 hours)
1. Create `src/algorithms/tinit.py`
2. Implement `TinitAlgorithm(BaseAlgorithm)`
3. Move algorithm-specific logic from `src/tinit.py`
4. Test: run both old and new, compare outputs

### Phase 4: Refactor CORDA (2 hours)
1. Create `src/algorithms/corda.py`
2. Implement `CordaAlgorithm(BaseAlgorithm)`
   - Note: CORDA uses different integration (no Troppo)
   - Override `load_expression()` for confidence mapping
3. Test: run both old and new, compare outputs

### Phase 5: Refactor Remaining Algorithms (2 hours)
1. FastCORE
2. RIPTiDe

### Phase 6: Update CLI/Flow (1 hour)
1. Update `scripts/flow.py` to use new classes
2. Update `scripts/cli.py` if needed
3. Test full pipeline end-to-end

---

## Testing Strategy

### Unit Tests
```python
# tests/test_algorithms.py
import pytest
from src.algorithms.gimme import GimmeAlgorithm
from src.algorithms.tinit import TinitAlgorithm

def test_base_algorithm_config_loading():
    """Test config loading in base class."""
    algo = GimmeAlgorithm()
    assert algo.cfg is not None

def test_model_type_detection():
    """Test model type auto-detection."""
    algo = GimmeAlgorithm()
    assert algo.model_type in ['human-gem', 'recon3d']

def test_expression_path_resolution():
    """Test expression file auto-detection."""
    algo = GimmeAlgorithm()
    assert algo.expr_path is not None
```

### Integration Tests
- Run old script → save output
- Run new class → save output
- Compare: reaction count, objective value, model size

---

## Benefits of This Refactoring

| Aspect | Before | After |
|--------|--------|-------|
| Config loading | Copy-pasted 5 times | Single method |
| Model path resolution | Copy-pasted 5 times | Single method |
| Expression auto-detection | Copy-pasted 5 times | Single method |
| Model export | Copy-pasted 5 times | Single method |
| New algorithm | ~150 lines boilerplate | ~30 lines |
| Bug fixes | Edit 5 files | Edit 1 file |

---

## Common Pitfalls to Avoid

1. **Don't over-abstract**: Only refactor what's actually duplicated
2. **Keep backward compatibility**: Ensure old CLI still works during transition
3. **Test incrementally**: Test each class after creation, not at the end
4. **Handle edge cases**: Some algorithms may need to override base methods (e.g., CORDA)

---

## File Naming Convention

After refactoring:
- `src/algorithms/` - All algorithm classes
- `src/scripts/` - Original scripts (keep for backward compatibility during transition)
- Move originals to `src/legacy/` after testing new classes

---

## Questions to Answer Before Starting

1. Should the original scripts be kept or replaced?
2. Should there be a factory pattern for algorithm instantiation?
3. How should algorithm-specific configuration be handled? (currently in `config.yaml`)
4. Should there be a unified CLI for all algorithms or keep separate entry points?

---

## Quick Reference Commands

```bash
# Create the directory structure
mkdir -p src/algorithms

# Run tests
pytest tests/test_algorithms.py -v

# Run a specific algorithm
python -c "from src.algorithms.gimme import GimmeAlgorithm; GimmeAlgorithm().run()"
```

---

## Key OOP Concepts Practiced

1. **Inheritance**: `GimmeAlgorithm(BaseAlgorithm)` inherits all common methods
2. **Abstraction**: `@abstractmethod` decorators enforce method implementation
3. **Encapsulation**: Data is stored in instance variables (`self.model`, `self.cfg`)
4. **Polymorphism**: Each algorithm implements `run()` differently
5. **Composition**: Uses external libraries (Cobra, Troppo) within class methods
6. **Factory Pattern** (optional future enhancement): `AlgorithmFactory.create('gimme')`
