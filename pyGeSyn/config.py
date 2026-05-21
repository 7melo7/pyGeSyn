import json
from pathlib import Path


def load_config(config_path):
    config = {}
    for cp in config_path.split(','):
        cp = cp.strip()
        with open(cp) as f:
            part = json.load(f)
        if not isinstance(part, dict):
            raise ValueError(f"Config in '{cp}' must be a JSON object")
        for name, paths in part.items():
            if name in config:
                existing = config[name]
                same = (existing.get('fasta') == paths.get('fasta') and
                        existing.get('gff3') == paths.get('gff3') and
                        existing.get('te') == paths.get('te'))
                if not same:
                    raise ValueError(
                        f"Duplicate genome '{name}' with different paths "
                        f"across config files")
                continue
            if 'fasta' not in paths:
                raise ValueError(f"Genome '{name}': missing 'fasta' key")
            p = Path(paths['fasta'])
            if not p.exists():
                raise FileNotFoundError(f"Genome '{name}': fasta not found: {p}")
            for key in ['gff3', 'te']:
                if key in paths:
                    p = Path(paths[key])
                    if not p.exists():
                        raise FileNotFoundError(f"Genome '{name}': {key} not found: {p}")
            config[name] = paths

    return config
