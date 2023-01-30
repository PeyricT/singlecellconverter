import argparse
from .converter import seurat2anndata, anndata2seurat, compress_h5ad


def _commands_seurat2anndata():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Seurat object path", required=True)
    parser.add_argument("-o", "--output", help="Anndata object path")
    parser.add_argument("-s", "--assay", help="Select which assay use for convertion default `auto`."
                                              " `all` will create an anndata files for each assay", default='auto')
    parser.add_argument("-g", "--graph", help="Select which graph use for convertion default `auto`.", default='auto')
    parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="store_true", default=False)
    parser.add_argument("-w", "--warnings", help="show warnings", action="store_true", default=False)
    args = parser.parse_args()

    if not args.input:
        raise ValueError(f'A Seurat object path is requiere [-i]')

    ext = args.input.split('.')[-1]
    args.input = args.input.replace(ext, ext.lower())
    if ext.lower() != 'h5seurat':
        raise ValueError(f'{args.input} is not a h5Seurat format file')

    seurat2anndata(args.input, args.output, args.assay, args.graph, args.verbosity, args.warnings)


def _commands_anndata2seurat():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Seurat object path", required=True)
    parser.add_argument("-o", "--output", help="Anndata object path")
    parser.add_argument("-n", "--nolayers", help="Only use X matrix, ignore layers", action="store_true", default=False)
    parser.add_argument("-c", "--counts", help="layers use for count, default 'X'", default='X')
    parser.add_argument("-d", "--data", help="layers use for data, default 'X'", default='X')
    parser.add_argument("-s", "--scale", help="layers use for scale.data, default None", default=None)
    parser.add_argument("-v", "--verbosity", help="Increase output verbosity", action="store_true", default=False)
    parser.add_argument("-w", "--warnings", help="Show warnings", action="store_true", default=False)
    args = parser.parse_args()

    if not args.input:
        raise ValueError(f'An Anndata object path is requiere [-i]')
    elif not args.input.endswith('.h5ad'):
        raise ValueError(f'{args.input} is not a h5ad format file')

    if not args.output:
        args.output = args.input.replace('.h5ad', '.h5Seurat')

    anndata2seurat(
        args.input,
        output_file=args.output,
        nolayers=args.nolayers,
        counts=args.counts,
        data=args.data,
        scale=args.scale,
        verbosity=args.verbosity,
        warnings=args.warnings,
    )


def _commands_compress_h5ad():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Seurat object path", required=True)
    parser.add_argument("-o", "--output", help="Anndata object path")
    parser.add_argument("-v", "--verbosity", help="Increase output verbosity", action="store_true", default=False)
    parser.add_argument("-w", "--warnings", help="Show warnings", action="store_true", default=False)
    args = parser.parse_args()

    if not args.input:
        raise ValueError(f'An Anndata object path is requiere [-i]')
    elif not args.input.endswith('.h5ad'):
        raise ValueError(f'{args.input} is not a h5ad format file')

    compress_h5ad(args.input, args.output, args.verbosity, args.warnings)
