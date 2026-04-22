"""
Immunex CLI Module

Command-line interface for Immunex tools.

Main entry point: imn (Immunex)
"""

from .main import main as imn_main


def batch_pdb_main(*args, **kwargs):
    """Lazy import batch_pdb entry point to avoid unrelated import failures."""
    from .batch_pdb import main
    return main(*args, **kwargs)


__all__ = ['imn_main', 'batch_pdb_main']
