"""Preserve intermediate files in a separate directory for each invocation."""
from pathlib import Path
import tempfile


def create_workdir():
    root = Path('temp')
    root.mkdir(exist_ok=True)
    return Path(tempfile.mkdtemp(prefix='pyGeSyn-', dir=root))
