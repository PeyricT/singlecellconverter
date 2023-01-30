from setuptools import setup, find_packages

setup(
    name='singlecellconverter',
    version='0.1.0',
    description='A tool for converting Anndata object to Seurat object and vice versa',
    url='https://github.com/PeyricT/singlecellconverter',
    author='Thibaut Peyric',
    author_email='thibaut.peyric@proton.me',
    license='BSD 2-clause',
    packages=find_packages(),
    install_requires=['h5py',
                      'numpy',
                      'pandas',
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.8',
    ],

    entry_points={
        'console_scripts': [
            'seurat2anndata = singlecellconverter._commands:_commands_seurat2anndata',
            'anndata2seurat = singlecellconverter._commands:_commands_anndata2seurat',
            'compressh5ad = singlecellconverter._commands:_commands_compress_h5ad',
        ],
    },
)
