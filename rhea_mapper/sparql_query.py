import csv
import os
import re
import json
import urllib.request

from rdflib import Graph
from SPARQLWrapper import SPARQLWrapper, JSON
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
	# Import to be compatible with biopython version lesser than 1.78
	from Bio.Alphabet.IUPAC import protein
except ImportError:
	# Exception to be compatible with biopython version superior to 1.78
	protein = None

from cobra import Model
from cobra.io import read_sbml_model, write_sbml_model

from rhea_mapper import rhea_reconstruction


def query_rhea(rhea_sparql_pathname, database_folder, output_folder):
	rhea_rdf_file = database_folder + '/rhea.rdf'
	g = Graph()
	g.parse(rhea_rdf_file)

	rhea_reg = r'[?a-zA-Z\_]*\srdfs:subClassOf\srh:Reaction'
	reg_expr = re.compile(rhea_reg)

	with open(rhea_sparql_pathname, 'r') as rhea_sparql_file:
		rhea_sparql_query = rhea_sparql_file.read()

	search_result = reg_expr.search(rhea_sparql_query)
	rhea_predicate = search_result.group().strip('').split(' ')[0].replace('?','')

	query_result = g.query(rhea_sparql_query)

	reactions = []
	for row in query_result:
		reactions.append(row[rhea_predicate].split('/')[-1])

	query_reactions = set(reactions)

	with open(output_folder+'/sparql_query_results.tsv', 'w') as tsv_file:
		csvwriter = csv.writer(tsv_file, delimiter='\t')
		csvwriter.writerow(['rhea_reaction'])
		for reaction_id in query_reactions:
			csvwriter.writerow([reaction_id])


def query_rhea_endpoint(rhea_sparql_pathname, output_folder):
	uniprot_sparql_endpoint = 'https://sparql.rhea-db.org/sparql'
	sparql = SPARQLWrapper(uniprot_sparql_endpoint)

	rhea_reg = r'[?a-zA-Z\_]*\srdfs:subClassOf\srh:Reaction'
	reg_expr = re.compile(rhea_reg)

	with open(rhea_sparql_pathname, 'r') as rhea_sparql_file:
		rhea_sparql_query = rhea_sparql_file.read()

	search_result = reg_expr.search(rhea_sparql_query)
	rhea_predicate = search_result.group().strip('').split(' ')[0].replace('?','')

	sparql.setQuery(rhea_sparql_query)

	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	query_reactions = []
	for result in results['results']['bindings']:
		protein_id = result[rhea_predicate]['value'].split('/')[-1]
		query_reactions.append(protein_id)

	with open(output_folder+'/sparql_query_results.tsv', 'w') as tsv_file:
		csvwriter = csv.writer(tsv_file, delimiter='\t')
		csvwriter.writerow(['rhea_reaction'])
		for reaction_id in query_reactions:
			csvwriter.writerow([reaction_id])


def rhea_sbml_creation(rhea_sbml_file, output_folder):
	rhea_model = read_sbml_model(rhea_sbml_file)

	query_reactions = []
	with open(output_folder+'/sparql_query_results.tsv', 'r') as tsv_file:
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
	write_sbml_model(species_model, output_folder+'/sparql_query.sbml')


def query_uniprot_protein(uniprot_sparql_query_pathanme, database_folder, output_folder, nb_cpu):
	uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
	sparql = SPARQLWrapper(uniprot_sparql_endpoint)

	uniprot_reg = r'[?a-zA-Z\_]*\sa\sup:Protein'
	reg_expr = re.compile(uniprot_reg)

	with open(uniprot_sparql_query_pathanme, 'r') as uniprot_sparql_file:
		uniprot_sparql_query = uniprot_sparql_file.read()

	search_result = reg_expr.search(uniprot_sparql_query)
	uniprot_predicate = search_result.group().strip('').split(' ')[0].replace('?','')

	sparql.setQuery(uniprot_sparql_query)

	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	proteins_ids = []
	for result in results['results']['bindings']:
		protein_id = result[uniprot_predicate]['value'].split('/')[-1]
		proteins_ids.append(protein_id)

	proteins_ids = set(proteins_ids)

	proteins_in_database = {record.id.split('|')[1]: record
							for record in SeqIO.parse(database_folder+'/uniprot_sprot.fasta', 'fasta')
							if record.id.split('|')[1] in proteins_ids}

	for record_id in proteins_in_database:
		proteins_in_database[record_id].id = proteins_in_database[record_id].id.split('|')[1]

	missing_proteins = list(set(proteins_ids) - set(list(proteins_in_database.keys())))

	if missing_proteins != []:
		uri_missing_proteins = ['up:'+prot_id for prot_id in missing_proteins]

		missing_proteins_query = """PREFIX up: <http://purl.uniprot.org/core/>
				PREFIX rh: <http://rdf.rhea-db.org/>
				PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

				SELECT distinct ?protein ?aminoacidsequence WHERE {{

					?protein owl:disjointWith ?sequence .
					?sequence rdfs:comment ?aminoacidsequence .
					# get EC associated to protein
					#?protein up:enzyme ?enzyme .
					# Get Rhea reaction linked to protein
					#?protein up:annotation ?a .
					#?a a up:Catalytic_Activity_Annotation .
					#?a up:catalyticActivity ?ca .
					#?ca up:catalyzedReaction ?reaction .
					VALUES ?protein {{ {0} }}
				}}""".format(' '.join(uri_missing_proteins))

		sparql.setQuery(missing_proteins_query)

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
	else:
		missing_proteins_records = []

	fasta_records = list(proteins_in_database.values()) + missing_proteins_records

	SeqIO.write(fasta_records, output_folder+'/sparql_query.fasta', 'fasta')

	rhea_reconstruction.manage_genome('sparql_query', None, output_folder+'/sparql_query.fasta', database_folder, output_folder, nb_cpu)
