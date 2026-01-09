import yaml

class Config:
    def __init__(self, config_dict):
        self._config = config_dict  # Store dict
    
    def __getattr__(self, key):
        return self._config.get(key)
    def __getitem__(self, key):
        return self._config[key]
    
    def get(self, key, default=None):
        return self._config.get(key, default) # with default value
    
    @classmethod
    def from_yaml(cls, path):
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(data)