import csv
import os
import re
import json
import urllib.request

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rdflib import Graph
from SPARQLWrapper import SPARQLWrapper, JSON

try:
	# Import to be compatible with biopython version lesser than 1.78
	from Bio.Alphabet.IUPAC import protein
except ImportError:
	# Exception to be compatible with biopython version superior to 1.78
	protein = None


from rhea_mapper import rhea_reconstruction


def query_sparql_uniprot_organism(organism, database_folder, output_folder, nb_cpu):
	organism_escaped = organism.replace(' ', '_')

	if not os.path.exists(output_folder+'/sparql_query/'):
		os.mkdir(output_folder+'/sparql_query/')

	if not os.path.exists(output_folder+'/annotation_evidence/'):
		os.mkdir(output_folder+'/annotation_evidence/')

	sparql_organism_query = output_folder+'/sparql_query/'+organism_escaped+'.rq'
	sparql_ec_query = output_folder+'/sparql_query/'+organism_escaped+'_ec_evidence.rq'
	sparql_rhea_query = output_folder+'/sparql_query/'+organism_escaped+'_rhea_evidence.rq'

	ec_evidence = output_folder+'/annotation_evidence/'+organism_escaped+'_ec_evidence.tsv'
	rhea_evidence = output_folder+'/annotation_evidence/'+organism_escaped+'_rhea_evidence.tsv'

	organism_ec = output_folder+'/'+organism_escaped+'_ec.tsv'
	organism_fasta = output_folder+'/'+organism_escaped+'.fasta'

	taxon_id = find_taxon_id(organism, database_folder)

	group_name_to_sparql_query(taxon_id, sparql_organism_query, sparql_ec_query, sparql_rhea_query)
	query_uniprot_annotation(sparql_ec_query, 'enzyme', ec_evidence)
	query_uniprot_annotation(sparql_rhea_query, 'reaction', rhea_evidence)

	prot_ecs = {}
	with open(ec_evidence, 'r') as tsv_file:
		csvreader = csv.reader(tsv_file, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			if line[0] not in prot_ecs:
				prot_ecs[line[0]] = [line[1]]
			else:
				prot_ecs[line[0]].append(line[1])

	with open(organism_ec, 'w') as tsv_file:
		csvwriter = csv.writer(tsv_file, delimiter='\t')
		csvwriter.writerow(['gene', 'EC_Number'])
		for prot_id in prot_ecs:
			csvwriter.writerow([prot_id, ','.join(prot_ecs[prot_id])])

	uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
	sparql = SPARQLWrapper(uniprot_sparql_endpoint)

	uniprot_reg = r'[?a-zA-Z\_]*\sa\sup:Protein'
	reg_expr = re.compile(uniprot_reg)

	with open(sparql_organism_query, 'r') as uniprot_sparql_file:
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

	fasta_records = list(proteins_in_database.values())

	SeqIO.write(fasta_records, organism_fasta, 'fasta')

	rhea_reconstruction.manage_genome(organism_escaped, organism_ec, organism_fasta, database_folder, output_folder, nb_cpu)


def find_taxon_id(input_taxon_name, database_folder):
	tax_id_matches = []

	taxon_names = {}
	unique_names = {}
	with open(database_folder+'/ncbi_taxonomy.dmp') as input_file:
		for line in input_file:
			tax_id, taxon_name, unique_name, _, _ = [element.strip('\t') for element in line.split('|')]
			taxon_names[taxon_name] = tax_id
			unique_names[unique_name] = tax_id

	if input_taxon_name in taxon_names:
		tax_id_matches.append(taxon_names[input_taxon_name])
	if input_taxon_name in unique_names:
		tax_id_matches.append(unique_names[input_taxon_name])

	if len(tax_id_matches) == 0:
		print('No taxon_id found for ' + input_taxon_name + ', look at ' + database_folder+'/ncbi_taxonomy.dmp' + ' to find the matching for your species.')

	return tax_id_matches[0]


def group_name_to_sparql_query(taxon_id, sparql_organism_query, sparql_ec_query, sparql_rhea_query):
	taxon_id_uri = 'taxon:' + taxon_id + ''

	organism_sparql_query = '''PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
	PREFIX up: <http://purl.uniprot.org/core/>

	SELECT DISTINCT ?protein
	{{
		# Taxonomy selection using the union of organism and taxonomy for the input taxon.
		{{
			?protein a up:Protein .
			?protein up:reviewed True .

			# Proteins from the targeted taxonomy/organism.
			?protein up:organism {0} .
		}}
		UNION
		{{
			?protein a up:Protein .
			?protein up:reviewed True .

			# Proteins from the targeted taxonomy/organism.
			?protein up:organism ?organism .
			?organism rdfs:subClassOf {0} .
		}}
	}}'''.format(taxon_id)

	ec_annotation_sparql_query = '''PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
	PREFIX up: <http://purl.uniprot.org/core/>

	SELECT DISTINCT ?protein ?enzyme ?evidence
	WHERE
	{{
		# Get for each protein the corresponding EC number with evidence.
		?protein up:enzyme ?enzyme .
		[] a rdf:Statement ;
					rdf:subject ?protein ;
					rdf:predicate up:enzyme ;
					rdf:object ?enzyme ;
					up:attribution ?attribution .
		?attribution up:evidence ?evidence .

		# Taxonomy selection using the union of organism and taxonomy for the input taxon.
		{{
			?protein a up:Protein .
			?protein up:reviewed True .

			# Proteins from the targeted taxonomy/organism.
			?protein up:organism {0} .
		}}
		UNION
		{{
			?protein a up:Protein .
			?protein up:reviewed True .

			# Proteins from the targeted taxonomy/organism.
			?protein up:organism ?organism .
			?organism rdfs:subClassOf {0} .
		}}
	}}'''.format(taxon_id)

	rhea_annotation_sparql_query = '''PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
	PREFIX up: <http://purl.uniprot.org/core/>

	SELECT DISTINCT ?protein ?reaction ?evidence
	{{
		# Get for each protein the corresponding Rhea reaction with evidence.
		?protein a up:Protein .
		?protein up:reviewed True ;
			up:annotation ?a ;
			up:attribution ?attribution .
		?a a up:Catalytic_Activity_Annotation ;
			up:catalyticActivity ?ca .
		?ca up:catalyzedReaction ?reaction .
		[] rdf:subject ?a ;
			rdf:predicate up:catalyticActivity ;
			rdf:object ?ca ;
			up:attribution ?attribution .
		?attribution up:evidence ?evidence .

		{{
			?protein a up:Protein .
			?protein up:reviewed True .

			# Proteins from the targeted taxonomy/organism.
			?protein up:organism {0} .
		}}
		UNION
		{{
			?protein a up:Protein .
			?protein up:reviewed True .

			# Proteins from the targeted taxonomy/organism.
			?protein up:organism ?organism .
			?organism rdfs:subClassOf {0} .
		}}
	}}'''.format(taxon_id)

	with open(sparql_organism_query, 'w') as output_file:
		output_file.write(organism_sparql_query)
	with open(sparql_ec_query, 'w') as output_file:
		output_file.write(ec_annotation_sparql_query)
	with open(sparql_rhea_query, 'w') as output_file:
		output_file.write(rhea_annotation_sparql_query)
	return True


def query_uniprot_annotation(annotation_query_pathname, query_annotation_type, output_file):
	uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
	sparql = SPARQLWrapper(uniprot_sparql_endpoint)

	with open(annotation_query_pathname, 'r') as annotation_query_file:
		annotation_query = annotation_query_file.read()
	sparql.setQuery(annotation_query)

	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	ec_evidence_results = []
	for result in results['results']['bindings']:
		protein_id = result['protein']['value'].split('/')[-1]
		protein_annotation = result[query_annotation_type]['value'].split('/')[-1]
		protein_annotation_evidence = result['evidence']['value'].split('/')[-1]

		ec_evidence_results.append([protein_id, protein_annotation, protein_annotation_evidence])

	with open(output_file, 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter='\t')
		csvwriter.writerow(['prot_id', query_annotation_type, 'evidence'])
		for prot_tuples in ec_evidence_results:
			prot_id = prot_tuples[0]
			ec_number = prot_tuples[1]
			evidence = prot_tuples[2]
			csvwriter.writerow([prot_id, ec_number, evidence])
