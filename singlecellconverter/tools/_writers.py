import h5py as hpy
import numpy as np


def write_matrix(m_adata, m_seurat, name, dim_=np.zeros(2)):
    if isinstance(m_seurat, hpy._hl.dataset.Dataset):
        m_adata.create_dataset(name, data=m_seurat, compression=True)

    elif isinstance(m_seurat, hpy._hl.group.Group):
        m_adata.create_group(name)
        m_adata[name].attrs['encoding-type'] = 'csr_matrix'
        m_adata[name].attrs['encoding-version'] = '0.1.0'
        m_adata[name].attrs['shape'] = dim_

        m_adata[name].create_dataset('data', data=m_seurat['data'], compression=True)
        m_adata[name].create_dataset('indices', data=m_seurat['indices'], compression=True)
        m_adata[name].create_dataset('indptr', data=m_seurat['indptr'], compression=True)

    else:
        raise TypeError(f"{type(m_seurat)} unknown, must be Group or Dataset")


def write_matrix_a2s(m_adata, m_seurat, name, dim_=np.zeros(2)):
    if isinstance(m_seurat, hpy._hl.dataset.Dataset):
        m_adata.create_dataset(name, data=m_seurat, compression=True)

    elif isinstance(m_seurat, hpy._hl.group.Group):
        m_adata.create_group(name)
        m_adata[name].attrs['dims'] = dim_

        m_adata[name].create_dataset('data', data=m_seurat['data'], compression=True)
        m_adata[name].create_dataset('indices', data=m_seurat['indices'], compression=True)
        m_adata[name].create_dataset('indptr', data=m_seurat['indptr'], compression=True)

    else:
        raise TypeError(f"{type(m_seurat)} unknown, must be Group or Dataset")


def write_metadata(m_adata, m_seurat, metadata_) -> None:
    if isinstance(m_seurat, hpy._hl.dataset.Dataset):
        m_adata.create_dataset(metadata_, data=m_seurat)

    elif isinstance(m_seurat, hpy._hl.group.Group):
        values = np.array(m_seurat['values'], dtype=np.int32) - 1
        m_adata.create_dataset(metadata_, data=values)
        m_adata['__categories'].create_dataset(metadata_, data=m_seurat['levels'])


def write_metadata_a2s(m_seurat, m_adata, metadata_) -> None:
    if '__categories' not in m_adata or metadata_ not in m_adata['__categories']:
        m_seurat.create_dataset(metadata_, data=m_adata[metadata_])
    else:
        m_seurat.create_group(metadata_)
        values = np.array(m_adata[metadata_], dtype=np.int32) + 1
        m_seurat[metadata_].create_dataset('levels', data=np.array(m_adata['__categories'][metadata_]))
        m_seurat[metadata_].create_dataset('values', data=values)


def write_reduction_dim(r_adata, r_seurat, reduction_):
    r_adata.create_dataset('X_' + reduction_, data=np.array(r_seurat['cell.embeddings']).transpose())


def write_reduction_dim_a2s(r_seurat, r_adata, reduction_):
    red = reduction_.replace('X_', '')
    r_seurat.create_group(red)
    r_seurat[red].create_dataset('cell.embeddings', data=np.array(r_adata['obsm'][reduction_]).transpose())
    if red == 'pca':
        features = r_adata['var']['_index']
        if 'higly_variable' in r_adata['var']:
            features = np.array(features)
            hvg = np.array(r_adata['var']['higly_variable'])
            features = features[hvg]

        r_seurat[red].create_dataset('features', data=features)


def write_group(g_adata, g_seurat, ignore=[]):
    for path in g_seurat:

        if path in ignore:
            continue

        if isinstance(g_seurat[path], hpy._hl.dataset.Dataset):
            comp = True
            data = g_seurat[path]
            if g_seurat[path].shape == () or g_seurat[path].shape is None:
                comp = None
            elif g_seurat[path].shape == (1,):
                comp = None
                data = np.array(data)[0]
            g_adata.create_dataset(path, data=data, compression=comp)
            for attr in g_seurat[path].attrs:
                if attr == 'categories':
                    continue
                g_adata[path].attrs[attr] = g_seurat[path].attrs[attr]

        elif isinstance(g_seurat[path], hpy._hl.group.Group):
            if path in g_adata:
                continue
            g_adata.create_group(path)
            for attr in g_seurat[path].attrs:
                g_adata[path].attrs[attr] = g_seurat[path].attrs[attr]
            write_group(g_adata[path], g_seurat[path])

        else:
            raise TypeError(f'`{type(g_seurat[path])}` is not handle.')


def write_misc_a2s(g_seurat, g_adata, ignore=[]):
    names = []
    for path in g_adata:
        if path in ignore:
            continue

        names.append(path)
        if isinstance(g_adata[path], hpy._hl.dataset.Dataset):
            comp = True
            data = g_adata[path]
            if g_adata[path].shape == () or g_adata[path].shape is None:
                comp = None
                data = [np.array(g_adata[path]).tolist()]
            g_seurat.create_dataset(path, data=data, compression=comp)
            for attr in g_adata[path].attrs:
                if attr == 'categories':
                    continue
                g_seurat[path].attrs[attr] = g_adata[path].attrs[attr]

        elif isinstance(g_adata[path], hpy._hl.group.Group):
            g_seurat.create_group(path)
            for attr in g_adata[path].attrs:
                g_seurat[path].attrs[attr] = g_adata[path].attrs[attr]
            write_misc_a2s(g_seurat[path], g_adata[path])

        else:
            raise TypeError(f'`{type(g_adata[path])}` is not handle.')

    g_seurat.attrs['names'] = names
