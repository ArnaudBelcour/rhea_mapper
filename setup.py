import os

from distutils.util import convert_path
from setuptools import setup, find_packages

init_data = {}
init_pathname = convert_path('rhea_mapper/__init__.py')
with open(init_pathname) as init_file:
    exec(init_file.read(), init_data)

setup(
    name='rhea_mapper',
    version=init_data['__version__'],
    packages=['rhea_mapper'],
    package_dir={'rhea_mapper': 'rhea_mapper'},
    entry_points={
        'console_scripts': [
            'rhea_mapper = rhea_mapper.__main__:main',
        ]
    },
    install_requires=['biopython', 'cobra', 'rdflib', 'sparqlwrapper'],
)