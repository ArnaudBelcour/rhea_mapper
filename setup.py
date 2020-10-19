from setuptools import setup

setup(
    name='rhea_mapper',
    packages=['rhea_mapper'],
    package_dir={'rhea_mapper': 'rhea_mapper'},
    entry_points={
        'console_scripts': [
            'rhea_mapper = rhea_mapper.__main__:main',
        ]
    },
    install_requires=['biopython', 'cobra', 'rdflib', 'sparqlwrapper'],
)