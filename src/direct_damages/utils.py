"""Utility functions."""

import yaml

def read_config(file):
    """Read the config yaml outside of snakemake workflow"""
    with open(file, 'r') as src:
        cfg = yaml.safe_load(src)
    return cfg