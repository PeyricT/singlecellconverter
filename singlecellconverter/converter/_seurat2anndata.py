import time as t
import h5py as hpy
import numpy as np
import pandas as pd
import warnings as ws

from .. import tools as tl


def seurat2anndata(input_file, output_file=None, assay_name='auto', graph_name='auto', verbosity=False, warnings=False):

    if warnings:
        ws.filterwarnings("default", category=Warning)
    else:
        ws.filterwarnings("ignore", category=Warning)

    if not output_file:
        output_file = input_file.replace('.h5seurat', '.h5ad')

    if verbosity:
        bar = tl.ProgressBar({
            'LOAD FILE                           ': 0,
            'TRANSFER MATRIX                     ': 6,
            'TRANSFER OBSERVATIONS               ': 50,
            'TRANSFER FEATURES                   ': 60,
            'TRANSFER DIMENSION REDUCTION MATRIX ': 70,
            'TRANSFER NEAREST NEIGHBOR GRAPHS    ': 80,
            'TRANSFER MARKER GENES               ': 90,
            'TRANSFER OTHERS                     ': 96,
            'DONE                                ': 100,
        })
        bar.start()

    # --- LOAD FILE --- #
    if verbosity:
        bar.update()

    seurat = hpy.File(input_file, 'r')

    if assay_name == 'auto':
        assay_name = next(iter(seurat['assays'].keys()))
    elif assay_name == 'all':
        assays = list(seurat['assays'].keys())
        seurat.close()
        for assay in assays:
            output_file_assay = output_file.replace('.h5ad', f'_{assay}.h5ad')
            seurat2anndata(input_file, output_file_assay, assay, verbosity, warnings)
        return

    adata = hpy.File(output_file, 'w')

    try:
        assay = seurat['assays'][assay_name]
    except KeyError:
        raise KeyError(f'`{assay_name}` assay is not in `{input_file}`')

    dim = tl.get_matrix_dim(assay['counts'])
    matrix = ['counts']

    # --- TRANSFER MATRIX --- #
    if verbosity:
        bar.update()

    if 'data' in assay and np.logical_and(*(tl.get_matrix_dim(assay['data']) == dim)):
        matrix.append('data')
    if 'scale.data' in assay and np.logical_and(*(tl.get_matrix_dim(assay['scale.data']) == dim)):
        matrix.append('scale.data')

    tl.write_matrix(adata, assay[matrix.pop()], 'X', dim)
    adata.create_group('layers')

    for mat in matrix:
        tl.write_matrix(adata['layers'], assay[mat], mat, dim)

    # --- TRANSFER OBSERVATIONS --- #
    if verbosity:
        bar.update()

    adata.create_group('obs')
    adata['obs'].create_group('__categories')

    columns = []
    for metadata in seurat['meta.data']:
        columns.append(metadata)
        tl.write_metadata(adata['obs'], seurat['meta.data'][metadata], metadata)

    adata['obs'].attrs['_index'] = '_index'
    adata['obs'].attrs['column-order'] = np.array(columns, dtype=object)
    adata['obs'].attrs['encoding-type'] = 'dataframe'
    adata['obs'].attrs['encoding-version'] = '0.1.0'

    # --- TRANSFER FEATURES --- #
    if verbosity:
        bar.update()

    adata.create_group('var')
    adata['var'].create_group('__categories')

    columns = []
    for metadata in assay['meta.features']:
        columns.append(metadata)
        tl.write_metadata(adata['var'], assay['meta.features'][metadata], metadata)

    adata['var'].attrs['_index'] = '_index'
    adata['var'].attrs['column-order'] = np.array(columns, dtype=object)
    adata['var'].attrs['encoding-type'] = 'dataframe'
    adata['var'].attrs['encoding-version'] = '0.1.0'

    if 'variable.features' in assay and 'highly_variable' not in adata['var']:
        hvg = np.array(assay['variable.features'])
        genes = pd.Series(assay['features'])

        adata['var'].create_dataset('highly_variable', data=genes.isin(hvg))
        adata['var'].attrs['column-order'] = np.concatenate((
            np.array(['highly_variable']),
            adata['var'].attrs['column-order']
        ))

    # ---  TRANSFER DIMENSION REDUCTION MATRIX --- #
    if verbosity:
        bar.update()

    adata.create_group('obsm')

    for reduction in seurat['reductions']:
        tl.write_reduction_dim(adata['obsm'], seurat['reductions'][reduction], reduction)

    # --- TRANSFER NEAREST NEIGHBOR GRAPHS --- #
    if verbosity:
        bar.update()

    adata.create_group('obsp')
    adata.create_group('uns')

    dim_obsp = np.array([dim[0], dim[0]])

    if graph_name == 'auto':
        graph_name = assay_name

    if graph_name+'_nn' in seurat['graphs'] and graph_name+'_snn' in seurat['graphs']:
        tl.write_matrix(adata['obsp'], seurat['graphs'][graph_name+'_nn'], 'distances', dim_obsp)
        tl.write_matrix(adata['obsp'], seurat['graphs'][graph_name+'_snn'], 'connectivities', dim_obsp)

        adata['uns'].create_group('neighbors')
        adata['uns']['neighbors'].create_dataset('distances_key', data='distances')
        adata['uns']['neighbors'].create_dataset('connectivities_key', data='connectivities')

        adata['uns']['neighbors'].create_group('params')
        adata['uns']['neighbors']['params'].create_dataset('method', data='seurat')
        adata['uns']['neighbors']['params'].create_dataset('metric', data='euclidean')
        adata['uns']['neighbors']['params'].create_dataset('n_neighbors', data='20')
        adata['uns']['neighbors']['params'].create_dataset('random_state', data=0)
    else:
        print(f'Warning : graphs {graph_name} not found in seurat object.')

    # --- TRANSFER MARKER GENES --- #
    if verbosity:
        bar.update()

    rank_test = list(seurat['misc']['marker_genes']['cerebro_seurat'].keys())
    if 'misc/adata/' in seurat:
        for rank in seurat['misc/adata']:
            if 'rank_'+rank in rank_test:
                rank_test.remove(rank)

    seurat_markers = seurat['misc']['marker_genes']['cerebro_seurat']
    for group in rank_test:
        adata['uns'].create_group('rank_' + group)

        cols = {}
        groupby = tl.get_categorical_series(seurat_markers[group][group])
        cols['logfoldchanges'] = np.array(seurat_markers[group]["avg_log2FC"])
        cols['names'] = np.array(seurat_markers[group]["gene"])
        cols['pvals'] = np.array(seurat_markers[group]["p_val"])
        cols['pvals_adj'] = np.array(seurat_markers[group]["p_val_adj"])

        max_ = groupby.value_counts()[0]

        df = pd.DataFrame([], columns=groupby.cat.categories)
        for col_name, col in cols.items():
            df = pd.DataFrame([], columns=groupby.cat.categories)
            for cat in groupby.cat.categories:
                tmp = col[groupby == cat]
                tmp2 = np.full((max_ - len(tmp)), np.nan)
                tmp = np.concatenate((tmp, tmp2))
                df[cat] = tmp

            data = df.to_records(index=False)
            if col_name == 'names':
                data = data.astype([(col, 'S16') for col in df.columns])
            adata['uns']['rank_' + group].create_dataset(col_name, data=data)

        df = -np.log10(df)
        df.replace([np.inf], 320, inplace=True)
        adata['uns']['rank_' + group].create_dataset('scores', data=df.to_records(index=False))

        adata['uns']['rank_' + group].create_group('params')
        adata['uns']['rank_' + group]['params'].create_dataset('groupby', data=group)
        adata['uns']['rank_' + group]['params'].create_dataset('reference', data="rest")
        adata['uns']['rank_' + group]['params'].create_dataset('use_raw', data=False)
        adata['uns']['rank_' + group]['params'].create_dataset('method', data='FindAllMarkers')

    # --- TRANSFER OTHERS --- #
    if verbosity:
        bar.update()

    adata['uns'].create_group('seurat')
    adata['uns']['seurat'].create_group('commands')
    adata['uns']['seurat'].create_group('misc')

    tl.write_group(adata['uns']['seurat']['commands'], seurat['commands'])
    tl.write_group(adata['uns']['seurat']['misc'], seurat['misc'], ignore=['adata'])
    if 'adata' in seurat['misc']:
        tl.write_group(adata['uns'], seurat['misc']['adata'])

    adata.close()
    seurat.close()

    if verbosity:
        bar.update()
        t.sleep(1)
