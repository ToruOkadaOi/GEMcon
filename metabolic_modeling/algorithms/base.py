from abc import ABC, abstractmethod
import cobra

# run, load and save #

class AlgorithmBase(ABC):
    def __init__(self, config):
        self.config = config
        self.model = None
        self.result = None
    
    @abstractmethod
    def run(self, model, **kwargs):
        pass
    
    def load_model(self, model_path: str):
        return cobra.io.read_sbml_model(model_path)
    
    def save_model(self, model, output_path: str):
        cobra.io.write_sbml_model(model, output_path)