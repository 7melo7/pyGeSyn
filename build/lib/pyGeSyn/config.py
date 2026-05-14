import json
from pathlib import Path


def load_config(config_path):
    with open(config_path) as f:
        config = json.load(f)

    if not isinstance(config, dict):
        raise ValueError("Config must be a JSON object mapping genome names to paths")

    for name, paths in config.items():
        for key in ['fasta', 'gff3', 'te']:
            if key not in paths:
                raise ValueError(f"Genome '{name}': missing '{key}' key")
            p = Path(paths[key])
            if not p.exists():
                raise FileNotFoundError(f"Genome '{name}': file not found: {p}")

    return config
