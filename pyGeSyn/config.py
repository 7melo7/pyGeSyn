import json
import re
from pathlib import Path


def load_gene_colors(path):
    """Read exact gene/transcript IDs paired with six-digit RGB colors."""
    colors = {}
    with open(path, encoding='utf-8-sig') as stream:
        for lineno, line in enumerate(stream, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = re.split(r'[,\s]+', line)
            if len(fields) != 2 or not re.fullmatch(r'#[0-9a-fA-F]{6}', fields[1]):
                raise ValueError(f'{path}:{lineno}: expected GeneID and #RRGGBB color')
            identifier, color = fields
            color = color.upper()
            if identifier in colors and colors[identifier] != color:
                raise ValueError(f'{path}:{lineno}: conflicting colors for {identifier}')
            colors[identifier] = color
    return colors


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
                        existing.get('te') == paths.get('te') and
                        existing.get('gene_colors') == paths.get('gene_colors'))
                if not same:
                    raise ValueError(
                        f"Duplicate genome '{name}' with different paths "
                        f"across config files")
                continue
            if 'fasta' not in paths:
                raise ValueError(f"Genome '{name}': missing 'fasta' key")
            if 'gene_colors' in paths and 'gff3' not in paths:
                raise ValueError(f"Genome '{name}': gene_colors requires gff3")
            p = Path(paths['fasta'])
            if not p.exists():
                raise FileNotFoundError(f"Genome '{name}': fasta not found: {p}")
            for key in ['gff3', 'te', 'gene_colors']:
                if key in paths:
                    p = Path(paths[key])
                    if not p.exists():
                        raise FileNotFoundError(f"Genome '{name}': {key} not found: {p}")
            config[name] = paths

    return config
