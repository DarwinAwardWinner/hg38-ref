'''Utility functions for use in Snakemake workflows.'''

import os
import os.path
import shutil

def ensure_empty_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
