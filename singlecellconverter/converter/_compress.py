import h5py as hpy
import numpy as np
import pandas as pd
import warnings as ws

from .. import tools as tl


def compress_h5ad(input_file, output_file=None, verbosity=False, warnings=False):

    if warnings:
        ws.filterwarnings("default", category=Warning)
    else:
        ws.filterwarnings("ignore", category=Warning)

    if not output_file:
        output_file = input_file.replace('.h5ad', '_compress.h5ad')

    default = hpy.File(input_file, 'r')
    comp = hpy.File(output_file, 'w')

    tl.write_group(comp, default)
    if verbosity:
        print('\rDone ...', sep='\r')

    comp.close()
    default.close()
