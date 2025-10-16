from setuptools import setup

setup(
    name='mapseq-processing',
    author='John Hover',
    version='0.9.0',
    include_package_data=True,
    description='MAPseq raw data processing.',
    packages=[
        'mapseq',
    ],
    package_data={
        '': ['*.hdf5']
    },
    install_requires=[
        'matplotlib',
        'scipy',
        'pandas',
    ]
)