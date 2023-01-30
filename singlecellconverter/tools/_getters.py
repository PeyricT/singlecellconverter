import h5py as hpy
import numpy as np
import pandas as pd


def get_matrix_dim(matrix_):
    if isinstance(matrix_, hpy._hl.dataset.Dataset):
        return matrix_.shape
    elif isinstance(matrix_, hpy._hl.group.Group):
        return np.flip(matrix_.attrs['dims'])
    else:
        raise TypeError(f"{type(matrix_)} unknown, must be Group or Dataset")


def get_matrix_dim_a2s(matrix_):
    if isinstance(matrix_, hpy._hl.dataset.Dataset):
        return np.flip(matrix_.shape)
    elif isinstance(matrix_, hpy._hl.group.Group):
        return np.flip(matrix_.attrs['shape'])
    else:
        raise TypeError(f"{type(matrix_)} unknown, must be Group or Dataset")


def get_categorical_array(group_):
    categories = {k: cat_.decode("utf-8") for k, cat_ in enumerate(group_['levels'])}
    values = group_['values']

    return np.array([categories[i - 1] for i in values])


def get_categorical_series(group_):
    categories = {k: cat_.decode("utf-8") for k, cat_ in enumerate(group_['levels'])}
    values = group_['values']

    return pd.Series([categories[i - 1] for i in values]).astype('category')