import csv
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cobra import Model
from cobra.io import read_sbml_model, write_sbml_model
from SPARQLWrapper import SPARQLWrapper, JSON

try:
    # Import to be compatible with biopython version lesser than 1.78
    from Bio.Alphabet.IUPAC import protein
except ImportError:
    # Exception to be compatible with biopython version superior to 1.78
    protein = None

from rhea_mapper import rhea_reconstruction

def rhea_sbml_from_file(element_file, rhea_sbml_file, output_folder):
	rhea_model = read_sbml_model(rhea_sbml_file)

	query_reactions = []
	with open(element_file, 'r') as tsv_file:
		csvreader = csv.reader(tsv_file, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			query_reactions.append(line[0])

	species_model = Model('test')

	sbml_reactions = []
	for reaction in rhea_model.reactions:
		if reaction.id.replace('R_','') in query_reactions:
			reaction.gene_reaction_rule = ''
			sbml_reactions.append(reaction)

	species_model.add_reactions(sbml_reactions)

	# Create sbml file.
	write_sbml_model(species_model, output_folder+'/list_rhea.sbml')

def uniprot_from_file(element_file, rhea_sbml_file, database_folder, output_folder, nb_cpu):
	proteins_ids = []
	with open(element_file, 'r') as tsv_file:
		csvreader = csv.reader(tsv_file, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			proteins_ids.append(line[0])

	proteins_ids = set(proteins_ids)

	proteins_in_database = {record.id: record
							for record in SeqIO.parse(database_folder+'/uniprot_sprot.fasta', 'fasta')
							if record.id in proteins_ids}

	missing_proteins = list(set(proteins_ids) - set(list(proteins_in_database.keys())))
	uri_missing_proteins = ['up:'+prot_id for prot_id in missing_proteins]

	uniprot_sparql_endpoint = 'https://sparql.rhea-db.org/sparql'
	sparql = SPARQLWrapper(uniprot_sparql_endpoint)

	sparql.setQuery("""PREFIX up: <http://purl.uniprot.org/core/>
            PREFIX rh: <http://rdf.rhea-db.org/>
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

            SELECT distinct ?protein ?aminoacidsequence WHERE {

				?protein owl:disjointWith ?sequence .
				?sequence rdfs:comment ?aminoacidsequence .
				VALUES ?protein { {0} }
            }""".format(' '.join(uri_missing_proteins)))

	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	missing_proteins_records = []
	for result in results['results']['bindings']:
		protein_id = result['protein']['value'].split('/')[-1]
		protein_seq = result['aminoacidsequence']['value']
		missing_proteins_records.append(Seq(protein_seq))
		if protein:
			fasta_record = SeqRecord(Seq(protein_seq, protein), id=protein_id, description='')
		else:
			fasta_record = SeqRecord(Seq(protein_seq), id=protein_id, description='')
		missing_proteins_records.append(fasta_record)

	fasta_records = list(proteins_in_database.values()) + missing_proteins_records

	SeqIO.write(fasta_records, output_folder+'/list_uniprot.fasta', 'fasta')

	rhea_reconstruction.manage_genome('list_uniprot', None, output_folder+'/list_uniprot.fasta', database_folder, output_folder, nb_cpu)
