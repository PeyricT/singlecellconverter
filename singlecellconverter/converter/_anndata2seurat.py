import h5py as hpy
import numpy as np
import pandas as pd
import warnings as ws

from .. import tools as tl


def anndata2seurat(
        input_file,
        output_file=None,
        nolayers=False,
        counts='X',
        data='X',
        scale=None,
        verbosity=False,
        warnings=False
):

    if warnings:
        ws.filterwarnings("default", category=Warning)
    else:
        ws.filterwarnings("ignore", category=Warning)

    if not output_file:
        output_file = input_file.replace('.h5ad', '.h5Seurat')

    if verbosity:
        bar = tl.ProgressBar({
            'LOAD FILE                  ': 0,
            'TRANSFER MATRIX            ': 6,
            'TRANSFER METAFEATURES      ': 50,
            'TRANSFER COMMANDS & MISC   ': 60,
            'TRANSFER GRAPHS            ': 70,
            'TRANSFER METADATA          ': 80,
            'TRANSFER REDUCTIONS        ': 90,
            'TRANSFER RANK TEST         ': 96,
            'DONE                       ': 100,
        })
        bar.start()

    # --- LOAD DATA --- #
    if verbosity:
        bar.update()

    adata = hpy.File(input_file, 'r')
    seurat = hpy.File(output_file, 'w')

    seurat.create_group('active.ident')
    seurat.create_group('assays')
    seurat.create_group('commands')
    seurat.create_group('graphs')
    seurat.create_group('meta.data')
    seurat.create_group('misc')
    seurat.create_group('reductions')
    seurat.create_group('images')
    seurat.create_group('tools')

    seurat.attrs['active.assay'] = ['RNA']
    seurat.attrs['project'] = ['scRNA']
    seurat.attrs['version'] = ['4.0.2']

    seurat['assays'].create_group('RNA')
    seurat['assays/RNA'].attrs['key'] = ['rna_']

    # --- WRITE MATRIX --- #
    if verbosity:
        bar.update()

    dim = tl.get_matrix_dim_a2s(adata['X'])

    for key, mat in dict(counts=counts, data=data).items():
        matrix = adata['layers'][mat] if mat in adata['layers'] else adata['X']
        tl.write_matrix_a2s(seurat['assays']['RNA'], matrix, key, dim)

    if scale and 'highly_variable' in adata['var']:
        matrix = np.array(adata['layers'][scale] if scale in adata['layers'] else adata['X'])
        matrix = matrix[np.array(adata['var']['highly_variable'])]

        tl.write_matrix_a2s(seurat['assays']['RNA'], matrix, 'scale.data')

    # --- WRITE ACTIVE INDENT --- #
    try:
        ident = next(iter(adata['obs']['__categories']))

        seurat['active.ident'].create_dataset('levels', data=adata['obs']['__categories'][ident])
        seurat['active.ident'].create_dataset('values', data=adata['obs'][ident])

    except StopIteration:
        pass

    # --- ASSAYS METAFEATURES --- #
    if verbosity:
        bar.update()

    seurat.create_dataset('cell.names', data=adata['obs']['_index'])
    seurat['assays']['RNA'].create_dataset('features', data=adata['var']['_index'])
    seurat['assays']['RNA'].create_group('meta.features')
    colnames_ = []

    if 'highly_variable' in adata['var']:
        hvg = np.array(adata['var']['highly_variable'])
        feat = np.array(adata['var']['_index'])
        seurat['assays']['RNA'].create_dataset('variable.features', data=feat[hvg])
        seurat['assays']['RNA']['meta.features'].create_dataset('vst.variable', data=hvg)
        colnames_.append('vst.variable')

        if scale:
            seurat['assays']['RNA'].create_dataset('scaled.features', data=feat[hvg])

    for feature in adata['var']:
        if feature == 'vst.variable' or feature == '__categories':
            continue
        colnames_.append(feature)
        tl.write_metadata_a2s(seurat['assays']['RNA']['meta.features'], adata['var'], feature)

    colnames_.remove('_index')
    seurat['assays']['RNA']['meta.features'].attrs['_index'] = ['_index']
    seurat['assays']['RNA']['meta.features'].attrs['colnames'] = colnames_

    if 'vst.variable' in colnames_:
        seurat['assays']['RNA']['meta.features'].attrs['logicals'] = ['vst.variable']

    # --- COMMANDS & MISC --- #
    if verbosity:
        bar.update()

    if 'seurat/commands' in adata['uns']:
        tl.write_misc_a2s(seurat['commands'], adata['uns']['seurat']['commands'])

    if 'seurat/misc' in adata['uns']:
        tl.write_misc_a2s(seurat['misc'], adata['uns']['seurat']['misc'])

    seurat['misc'].create_group('adata')
    tl.write_misc_a2s(seurat['misc']['adata'], adata['uns'], ignore=['seurat'])

    # --- GRAPHS --- #
    if verbosity:
        bar.update()

    dim_graph = np.array([dim[1], dim[1]])
    if 'connectivities' in adata['obsp']:
        tl.write_matrix_a2s(seurat['graphs'], adata['obsp']['connectivities'], 'RNA_snn', dim_graph)
    if 'distances' in adata['obsp']:
        tl.write_matrix_a2s(seurat['graphs'], adata['obsp']['distances'], 'RNA_nn', dim_graph)

    seurat['/graphs/RNA_nn'].attrs['assay.used'] = ['RNA']
    seurat['/graphs/RNA_snn'].attrs['assay.used'] = ['RNA']

    seurat['/graphs/RNA_nn'].attrs['dims'] = dim_graph
    seurat['/graphs/RNA_snn'].attrs['dims'] = dim_graph

    # --- METADATA --- #
    if verbosity:
        bar.update()

    colnames_ = []
    for metadata in adata['obs']:
        if metadata == '__categories':
            continue
        colnames_.append(metadata)
        tl.write_metadata_a2s(seurat['meta.data'], adata['obs'], metadata)

    seurat['meta.data'].attrs['_index'] = ['_index']
    seurat['meta.data'].attrs['colnames'] = colnames_

    # --- REDUCTIONS --- #
    if verbosity:
        bar.update()

    for reduction in adata['obsm']:
        reduc = reduction.replace('X_', '')
        tl.write_reduction_dim_a2s(seurat['reductions'], adata, reduction)
        seurat['reductions'][reduc].attrs['active.assay'] = ['RNA']
        if reduction == 'X_umap':
            seurat['reductions'][reduc].attrs['global'] = np.array([1])
            seurat['reductions'][reduc].attrs['key'] = ['UMAP_']
        elif reduction == 'X_pca':
            seurat['reductions'][reduc].attrs['global'] = np.array([0])
            seurat['reductions'][reduc].attrs['key'] = ['PC_']
        elif reduction == 'X_UMAP_3D':
            seurat['reductions'][reduc].attrs['global'] = np.array([1])
            seurat['reductions'][reduc].attrs['key'] = ['UMAP3D_']

    # --- RANK TEST --- #
    if verbosity:
        bar.update()

    rank_test = [
        np.array(adata[f'uns/{rank}/params/groupby']).astype(str).tolist() for rank in adata['uns'] if rank.startswith('rank_')
    ]
    path_test = [
        rank for rank in adata['uns'] if rank.startswith('rank_')
    ]

    if 'uns/seurat/misc/marker_genes/cerebro_seurat' in adata:
        for rank in adata['uns/seurat/misc/marker_genes/cerebro_seurat']:
            if rank in rank_test:
                rank_test.remove(rank)

    if 'marker_genes' not in seurat['misc']:
        seurat['misc'].create_group('marker_genes')
    if 'cerebro_seurat' not in seurat['misc/marker_genes']:
        seurat['misc/marker_genes'].create_group('cerebro_seurat')

    colnames_ = dict(
        logfoldchanges='avg_log2FC',
        pvals='p_val',
        pvals_adj='p_val_adj',
        scores='score',
        names='gene'
    )
    for rank, path in zip(rank_test, path_test):
        seurat['misc/marker_genes/cerebro_seurat'].create_group(rank)

        for adata_col, seurat_col in colnames_.items():
            col = np.array(adata['uns'][path][adata_col])
            col = col[:250]
            col = np.array([col[type_] for type_ in col.dtype.names])
            col = col.flatten()
            seurat[f'misc/marker_genes/cerebro_seurat/{rank}'].create_dataset(seurat_col, data=col)

        slice = []
        for k in range(len(np.array(adata['obs']['__categories'][rank]))):
            slice += [k+1]*250

        seurat[f'misc/marker_genes/cerebro_seurat/{rank}'].create_group(rank)
        seurat[f'misc/marker_genes/cerebro_seurat/{rank}/{rank}'].create_dataset('levels',
                                                                                 data=adata['obs']['__categories'][
                                                                                     rank])
        seurat[f'misc/marker_genes/cerebro_seurat/{rank}/{rank}'].create_dataset('values',
                                                                                 data=slice)
    if verbosity:
        bar.update()

    seurat.close()
    adata.close()
            


