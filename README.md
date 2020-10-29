# rhea_mapper

Rheamapper is a Python3 (Python >= 3.6) tool to create draft metabolic network from genomes, sparql query on Rhea or Uniprot and from tsv files containing list of Rhea reactions or Uniprot protein IDs.

## Table of contents
- [rhea_mapper](#rhea_mapper)
  - [Table of contents](#table-of-contents)
  - [License](#license)
  - [Requirements](#requirements)
  - [Installation with pip](#installation-with-pip)
  - [Usage](#usage)
    - [database](#database)
    - [genome](#genome)
    - [sparql](#sparql)
    - [taxonomy](#taxonomy)
    - [from_file](#from_file)

## License

This project is licensed under the GNU General Public License - see the [LICENSE](https://github.com/ArnaudBelcour/rhea_mapper/blob/master/LICENSE) file for details.


## Requirements

Python 3 (Python 3.6 is tested). rhea_mapper uses a certain number of Python dependencies:

* [biopython](https://github.com/biopython/biopython) to handle fasta files.
* [cobra](https://github.com/opencobra/cobrapy) to create sbml files.
* [rdflib](https://github.com/RDFLib/rdflib) to load rhea.rdf and query it.
* [sparqlwrapper](https://github.com/RDFLib/sparqlwrapper) to query Uniprot and Rhea sparql endpoint.

And another tool:
* [OrthoFinder](https://github.com/davidemms/OrthoFinder) to find orthologs to uniprot with experimental  evidence linked to rhea. 


## Installation with pip

```
pip install rhea_mapper
```

## Usage

rhea_mapper has 4 main commands: genome, sparql, database and from_file.

````
usage: rhea_mapper [-h] [-v] {database,genome,sparql,taxonomy,from_file} ...

Create Rhea sbml files from genomes. For specific help on each subcommand use:
rhea_mapper {cmd} --help

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  valid subcommands:

  {database,genome,sparql,taxonomy,from_file}
    database            Create rhea_mapper database [internet connection
                        required]
    genome              metabolic network reconstruction from genome
    sparql              Using a SPARQL query for reaction ID on Rhea or
                        protein ID on Uniprot, create the corresponding
                        metabolic network
    taxonomy            With a taxonomy, query Uniprot to create the
                        corresponding metabolic network
    from_file           Using a tsv file with list of Rhea reaction or Uniprot
                        protein, create the corresponding metabolic network

Requires: Orthofinder for genome if you provide a fasta
````

### database

This command downloads and creates the files needed by rhea_mapper to work. As it downlaods files, this funciton requires an internet connection.

By using the RDF file of Rhea, rhea_mapper will create a SBMl file, with some specificity:

- the direction of reactions can be: arbitrary, follow MetaCyc reaction direction or follow KEGG reaction direction.
- for the variable stoichiometric coefficient (n, n-1, n+1, 2n) in reaction, n has been set to 2.
- genes linked to Rhea reactions with experimental evidence are associated to these reactions.

````
usage: rhea_mapper database [-h] -d DATABASE

Create rhea_mapper database [internet connection required]

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Database folder, if not existing it will be created
````

It will create an output_folder containing:

- rhea.rdf: the RDF file of Rhea database.
- rhea2ec.tsv: mapping between Rhea reaction and EC number.
- rhea2uniprot_sprot.tsv: mapping between Rhea reaction and Uniprot protein ID.
- uniprot_sprot.fasta: fasta containg all amino-acids sequences from Swissprot.
- uniprot_rhea_evidence.tsv: tsv file containing Uniprot protein ID linked to Rhea reaction with experimental evidence.
- uniprot_rhea_evidence.fasta: fasta file containing protein linked to Rhea reaction with experimental evidence.
- rhea.sbml: SBML file created by rhea_mapper using the rhea.rdf file and the uniprot_rhea_evidence.tsv (to have gene association).

Rhea files come from the [Rhea download page](ftp://ftp.expasy.org/databases/rhea/rdf/).

uniprot_sprot.fasta is downloaded from the [Uniprot download page](https://www.uniprot.org/downloads). It is the `Reviewed (Swiss-Prot)`.

uniprot_rhea_evidence.* files are created with a SPARQL query on the [Uniprot SPARQL endpoint](https://sparql.uniprot.org/sparql).

### genome

This genome takes as input folders with subfolders containing a tsv file with annotation (gene linked to EC number) and a fasta file (proteome of the species of interest).

Then it will map the EC to Rhea reactions.

in a second step it will take the proteome to search for orthologs against a database of Uniprot protein linked to Rhea reactions with experimental evidence.

Then a sbml file will be created.

````
usage: rhea_mapper genome [-h] [-f FOLDER] -o OUTPUT -d DATABASE [-c CPU]

metabolic network reconstruction from genome

optional arguments:
  -h, --help            show this help message and exit
  -f FOLDER, --folder FOLDER
                        Input folder containing either annotation, genomes or
                        both
  -o OUTPUT, --output OUTPUT
                        Output folder
  -d DATABASE, --database DATABASE
                        Database folder, if not existing it will be created
  -c CPU, --cpu CPU     cpu number for multi-process
````

The structure of the input_folder is:

````
    input_folder
        ├── organism_1
        │   └── organism_1.gbk
        ├── organism_2
        │   └── organism_2.tsv
        │   └── organism_2.fasta
        ├── organism_3
        │   └── organism_3.tsv
        ├── organism_4
        │   └── organism_4.fasta
        ..
        └── organism_n
            └── organism_n.gbk
````

GenBank files will be converted into the format use by `rhea_mapper genome` meaning a tsv file + a fasta file. It is also possible to give only a fasta file or a tsv file.

The fasta file must contain amino-acids sequences.

The tsv file must be structured like:

| gene   | Ec_Number |
|--------|-----------|
| gene_1 | X.X.X.X   |
| ...    |           |
| gene_n | X.X.X.X   |

### sparql

Using a sparql query, rhea_mapper will query either a rhea rdf file from its database, the uniprot sparql endpoint or the rhea sparql endpoint.

A Uniprot query must contain the predicate of protein Id from uniprot (`?predicate a up:Protein`). The set of proteins find by the sparql query will then be used to create a draft metabolic network. Rhea_mapper will use the reconstruction method of `genome` to have a sbml file.

For the Rhea query, it must a set of Rhea reactions (`?predicate rdfs:subClassOf rh:Reaction`). Then using these IDs, rhea_mapper will create the corresponding sbml files.

````
usage: rhea_mapper sparql [-h] [-s SPARQL] -e ENDPOINT -d DATABASE -o OUTPUT
                          [-c CPU]

Using a SPARQL query for reaction ID on Rhea or protein ID on Uniprot, create
the corresponding metabolic network

optional arguments:
  -h, --help            show this help message and exit
  -s SPARQL, --sparql SPARQL
                        SPARQL query file
  -e ENDPOINT, --endpoint ENDPOINT
                        Endpoint either 'uniprot' or 'rhea' or 'rhea_endpoint
  -d DATABASE, --database DATABASE
                        Database folder, if not existing it will be created
  -o OUTPUT, --output OUTPUT
                        Output folder
  -c CPU, --cpu CPU     cpu number for multi-process
````

### taxonomy

From an organism name or a taxonomy group, rhea_mapper will query Uniprot to retrieve the proteins associated to this organism/taxonomy and will create a draft metabolic network form it.

````
usage: rhea_mapper taxonomy [-h] -d DATABASE -o OUTPUT [-c CPU]
                            [--organism ORGANISM] [--taxonomy TAXONOMY]

With a taxonomy, query Uniprot to create the corresponding metabolic network

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Database folder, if not existing it will be created
  -o OUTPUT, --output OUTPUT
                        Output folder
  -c CPU, --cpu CPU     cpu number for multi-process
  --organism ORGANISM   Name of organism or group of organism to query uniprot
                        for Protein
  --taxonomy TAXONOMY   TSV file with taxonomy name of organism
````

### from_file

Using a tsv file with either Rhea reactions ID or Uniprot protein IDs, rhea_mapper will create draft metabolic network.

For the Rhea reactions, it will create a sbml file.

For the Uniprot protein, it will use the reconstruction method of `genome` to have a sbml file.

````
usage: rhea_mapper from_file [-h] -l LIST -d DATABASE -o OUTPUT [-c CPU] -r
                             REFERENCE

Using a tsv file with list of Rhea reaction or Uniprot protein, create the
corresponding metabolic network

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --l LIST     File containing list of Rhea reacitons IDs or protein
                        IDs
  -d DATABASE, --database DATABASE
                        Database folder, if not existing it will be created
  -o OUTPUT, --output OUTPUT
                        Output folder
  -c CPU, --cpu CPU     cpu number for multi-process
  -r REFERENCE, --reference REFERENCE
                        Database for reference either 'uniprot' or 'rhea'
````
