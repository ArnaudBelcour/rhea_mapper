import csv
import gzip
import json
import os
import pkg_resources
import sys
import shutil
import urllib.request
import zipfile

from Bio import SeqIO
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model
from datetime import datetime
from rdflib import Graph
from SPARQLWrapper import SPARQLWrapper, JSON

from Bio import __version__ as biopython_version
from cobra import __version__ as cobra_version
from rdflib import __version__ as rdflib_version
from SPARQLWrapper import __version__ as sparqlwrapper_version


def urllib_reporthook(count, block_size, total_size):
	downloaded_size = int(count * block_size) / (1024 * 1024)
	percent = min(int(count*block_size*100/total_size),100)
	total_size_show = total_size / (1024 * 1024)
	sys.stdout.write("\rDownloading: {0}%, {1:.2f} MB on {2:.2f} MB".format(percent, downloaded_size, total_size_show))
	sys.stdout.flush()


def rhea_to_sbml(rhea_rdf_file, uniprot_rhea_evidence, output_rhea_sbml, output_rhea_variable_stoichiometric_reactions):
	# Read rhea rdf file.
	g = Graph()
	g.parse(rhea_rdf_file)

	# Get all the approved reactions with their metabolites.
	query_reactions_metabolites = g.query(
		"""PREFIX rh:<http://rdf.rhea-db.org/>
		PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
		SELECT ?reaction ?reactionSide ?compound ?compoundName ?compoundFormula ?coefficient
		WHERE {{
		?reaction rdfs:subClassOf rh:Reaction .
		?reaction rh:status rh:Approved .
		?reaction rh:side ?reactionSide .
		?reactionSide rh:contains ?participant .
		?reactionSide ?rhContains ?participant .
		?rhContains rdfs:subPropertyOf rh:contains .
		?rhContains rh:coefficient ?coefficient .
		?participant rh:compound ?compound .
		?compound rh:name ?compoundName .
		OPTIONAL
			{{
			?compound rh:formula ?compoundFormula .
			}}
		}}""")

	# Find Uniprot protein associated to Rhea reactions.
	rhea_uniprots = {}
	with open(uniprot_rhea_evidence, 'r') as rhea_file:
		csvreader = csv.reader(rhea_file, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			reactions = line[1].split(',')
			protein_id = line[0]
			for reaction in reactions:
				if reaction not in rhea_uniprots:
					rhea_uniprots[reaction] = [protein_id]
				else:
					rhea_uniprots[reaction].append(protein_id)

	# Get all xrefs linked to reversible reactions.
	query_reversible_reaction = g.query(
		"""PREFIX rh:<http://rdf.rhea-db.org/>
		PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
		SELECT ?reaction ?bidirection ?xref
		WHERE {{
		?reaction rdfs:subClassOf rh:Reaction .
		?reaction rh:status rh:Approved .
		?reaction rh:side ?reactionSide .
		?reaction rh:bidirectionalReaction ?bidirection .
		?bidirection rdfs:seeAlso ?xref .
		}}""")

	# Extract xref linked to bidrections.
	# This way we have all reactions with reactions from other database indicating that they are reversible.
	reversible_reactions = {}
	for row in query_reversible_reaction:
		reaction_id, bireaction_id, bireaction_xref = row
		reaction_id_sbml = reaction_id.split('/')[-1]
		if reaction_id_sbml not in reversible_reactions:
			reversible_reactions[reaction_id_sbml] = [bireaction_xref]
		else:
			reversible_reactions[reaction_id_sbml].append(bireaction_xref)

	# Get all reaction with one direction.
	query_directed_reaction = g.query(
		"""PREFIX rh:<http://rdf.rhea-db.org/>
		PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
		SELECT ?reaction ?compoundproduct ?compoundsubstrate ?xref
		WHERE {{
		?reaction rdfs:subClassOf rh:Reaction .
		?reaction rh:status rh:Approved .
		?reaction rh:side ?reactionSide .
		?reaction rh:directionalReaction ?direaction .
		?direaction rh:products ?product .
		?direaction rh:substrates ?substrate .
		?product rh:contains ?participantproduct .
		?participantproduct rh:compound ?compoundproduct .
		?substrate rh:contains ?participantsubstrate .
		?participantsubstrate rh:compound ?compoundsubstrate .
		?direaction rdfs:seeAlso ?xref .
		}}""")

	# Get all the reactions with direction except the reversible ones.
	directed_reactions = {}
	for row in query_directed_reaction:
		reaction_id, reaction_product, reaction_substrate, reaction_xref = row
		reaction_product = reaction_product.split('/')[-1]
		reaction_substrate = reaction_substrate.split('/')[-1]
		reaction_id_sbml = reaction_id.split('/')[-1]
		if reaction_id_sbml not in reversible_reactions:
			if reaction_id_sbml not in directed_reactions:
				directed_reactions[reaction_id_sbml] = {}
				directed_reactions[reaction_id_sbml]["xrefs"] = [reaction_xref]
				directed_reactions[reaction_id_sbml]["product"] = [reaction_product]
				directed_reactions[reaction_id_sbml]["substrate"] = [reaction_substrate]
			else:
				if reaction_xref not in directed_reactions[reaction_id_sbml]["xrefs"]:
					directed_reactions[reaction_id_sbml]["xrefs"].append(reaction_xref)
				if reaction_product not in directed_reactions[reaction_id_sbml]["product"]:
					directed_reactions[reaction_id_sbml]["product"].append(reaction_product)
				if reaction_substrate not in directed_reactions[reaction_id_sbml]["substrate"]:
					directed_reactions[reaction_id_sbml]["substrate"].append(reaction_substrate)

	# Create cobra model.
	model = Model('rhea')

	# Create cobra metabolites dictionary for each reactions.
	reactions = {}
	metabolites = {}
	variable_stoich_reactions = []
	for row in query_reactions_metabolites:
		reaction_id, reaction_side, metabolite_id, metabolite_name, metabolite_formula, stoichiometric_coefficient = row
		if reaction_id not in reactions:
			reaction_id_sbml = reaction_id.split('/')[-1]
			metabolite_id = metabolite_id.split('/')[-1].replace('Compound_', '')
			# Create metabolite.
			if metabolite_id not in metabolites:
				metabolite_id_sbml = Metabolite(metabolite_id, compartment='c', name=metabolite_name, formula=metabolite_formula)
				metabolites[metabolite_id] = metabolite_id_sbml
			else:
				metabolite_id_sbml = metabolites[metabolite_id]
			if reaction_id_sbml not in reactions:
				reactions[reaction_id_sbml] = {}

			# Manage specific stoechiometric coefficient.
			if str(stoichiometric_coefficient) in ['N', 'Nplus1', 'Nminus1', '2n']:
				variable_stoich_reactions.append([reaction_id_sbml, metabolite_id, str(stoichiometric_coefficient)])
			if str(stoichiometric_coefficient) == 'N':
				stoichiometric_coefficient = 2.0
			elif str(stoichiometric_coefficient) == 'Nplus1':
				stoichiometric_coefficient = 3.0
			elif str(stoichiometric_coefficient) == 'Nminus1':
				stoichiometric_coefficient = 1.0
			elif str(stoichiometric_coefficient) == '2n':
				stoichiometric_coefficient = 4.0
			else:
				stoichiometric_coefficient = float(stoichiometric_coefficient)

			# Create reation in the direction find by the directed reaction SPARQL query.
			if reaction_id_sbml in directed_reactions:
				# Product compounds are assigned to 1.0 in cobra model.
				if metabolite_id in directed_reactions[reaction_id_sbml]['product']:
					reactions[reaction_id_sbml][metabolite_id_sbml] = stoichiometric_coefficient
				# Substrate compounds are assigned to -1.0 in cobra model.
				elif metabolite_id in directed_reactions[reaction_id_sbml]['substrate']:
					reactions[reaction_id_sbml][metabolite_id_sbml] = - stoichiometric_coefficient
			#Otherwise use the arbitrary direction.
			else:
				# Right compounds are assigned to 1.0 in cobra model.
				if reaction_side.endswith('_R'):
					reactions[reaction_id_sbml][metabolite_id_sbml] = stoichiometric_coefficient
				# Left compounds are assigned to -1.0 in cobra model.
				elif reaction_side.endswith('_L'):
					reactions[reaction_id_sbml][metabolite_id_sbml] = - stoichiometric_coefficient

	if len(variable_stoich_reactions) > 0:
		var_rxn_ids = []
		with open(output_rhea_variable_stoichiometric_reactions, 'w') as variable_stoich_file:
			csvwriter = csv.writer(variable_stoich_file, delimiter='\t')
			csvwriter.writerow(['reaction_id', 'metabolite_id', 'stoichiometric_coefficient'])
			for list_var_stoch_reaction in variable_stoich_reactions:
				csvwriter.writerow(list_var_stoch_reaction)
				var_rxn_ids.append(list_var_stoch_reaction[0])
		print('Warning there is ' + str(len(set(var_rxn_ids))) + ' reactions with variable stoichiometric coefficient.')
		print('The list of these reactions can be found in: ' + output_rhea_variable_stoichiometric_reactions)

	# Create each reaction and add their compounds.
	sbml_reactions = []
	for reaction_id in reactions:
		if reaction_id in reversible_reactions:
			# Set reaction reversible.
			reaction = Reaction(reaction_id, upper_bound=1000, lower_bound=-1000)
		else:
			reaction = Reaction(reaction_id)
		if reaction_id in rhea_uniprots:
			reaction.gene_reaction_rule = '( ' + ' or '.join([prot for prot in rhea_uniprots[reaction_id]]) + ' )'
		reaction.add_metabolites(reactions[reaction_id])

		sbml_reactions.append(reaction)

	model.add_reactions(sbml_reactions)

	# Create sbml file.
	write_sbml_model(model, output_rhea_sbml)


def download_database(database_folder):
	if not os.path.exists(database_folder):
		os.mkdir(database_folder)

	# Get Rhea release version
	rhearesponse = urllib.request.urlopen('ftp://ftp.expasy.org/databases/rhea/rhea-release.properties')
	rhea_lines = rhearesponse.readlines()
	rhea_release_number = rhea_lines[0].decode('utf-8').split('=')[1].replace('\n','')
	rhea_release_date = rhea_lines[1].decode('utf-8').split('=')[1].replace('\n','')

	# Get Uniprot release version
	uniprot_response = urllib.request.urlopen('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt')
	uniprot_lines = uniprot_response.readlines()
	swissprot_release_number = uniprot_lines[1].decode('utf-8').split(' ')[2].replace('\n','')
	swissprot_release_date = uniprot_lines[1].decode('utf-8').split(' ')[4].replace('\n','')
	trembl_release_number = uniprot_lines[2].decode('utf-8').split(' ')[2].replace('\n','')
	trembl_release_date = uniprot_lines[2].decode('utf-8').split(' ')[4].replace('\n','')

	print('Download Rhea RDF file')
	urllib.request.urlretrieve('ftp://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz', database_folder + '/rhea.rdf.gz', reporthook=urllib_reporthook)
	with gzip.open(database_folder + '/rhea.rdf.gz', 'rb') as f_in:
		with open(database_folder + '/rhea.rdf', 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	os.remove(database_folder + '/rhea.rdf.gz')
	print('\n')

	print('Download Rhea rhea2ec file')
	urllib.request.urlretrieve('ftp://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv', database_folder + '/rhea2ec.tsv', reporthook=urllib_reporthook)
	print('\n')

	print('Download Rhea rhea2uniprot file')
	urllib.request.urlretrieve('ftp://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot%5Fsprot.tsv', database_folder + '/rhea2uniprot_sprot.tsv', reporthook=urllib_reporthook)
	print('\n')

	print('Download Reviewed (Swiss-Prot) fasta file')
	urllib.request.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz', database_folder + '/uniprot_sprot.fasta.gz', reporthook=urllib_reporthook)
	with gzip.open(database_folder + '/uniprot_sprot.fasta.gz', 'rb') as f_in:
		with open(database_folder + '/uniprot_sprot.fasta', 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	os.remove(database_folder + '/uniprot_sprot.fasta.gz')
	print('\n')

	print('Download Reviewed NCBI taxonomy file')
	urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip', database_folder + '/taxdmp.zip', reporthook=urllib_reporthook)

	with zipfile.ZipFile(database_folder+'/taxdmp.zip',"r") as zip_taxonomy:
		zip_taxonomy.extract('names.dmp', database_folder)
		os.rename(database_folder+'/names.dmp', database_folder+'/ncbi_taxonomy.dmp')
	os.remove(database_folder + '/taxdmp.zip')
	print('\n')

	sparql = SPARQLWrapper('https://sparql.uniprot.org/sparql')

	print('Find Uniprot protein with evidence for Rhea reaction')
	sparql.setQuery("""PREFIX up: <http://purl.uniprot.org/core/>
			PREFIX rh: <http://rdf.rhea-db.org/>
			PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

			SELECT DISTINCT ?rhea_reaction ?protein WHERE {{

			# ECO 269 is experimental evidence
			BIND (<http://purl.obolibrary.org/obo/ECO_0000269> as ?evidence)
			?protein up:reviewed true ;
				up:annotation ?a ;
				up:attribution ?attribution .
			?a a up:Catalytic_Activity_Annotation ;
				up:catalyticActivity ?ca .
			?ca up:catalyzedReaction ?rhea_reaction .

			[] rdf:subject ?a ;
				rdf:predicate up:catalyticActivity ;
				rdf:object ?ca ;
				up:attribution ?attribution .
			?attribution up:evidence ?evidence .
			}}""")

	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	# Create mapping file uniprot ID -> Rhea Id with evidence.
	mapping_prot_rhea = {}
	for result in results["results"]["bindings"]:
		protein_id = result["protein"]["value"].split('/')[-1]
		reaction_id = result["rhea_reaction"]["value"].split('/')[-1]
		if protein_id not in mapping_prot_rhea:
			mapping_prot_rhea[protein_id] = [reaction_id]
		else:
			mapping_prot_rhea[protein_id].append(reaction_id)

	with open(database_folder + '/uniprot_rhea_evidence.tsv', 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter ='\t')
		csvwriter.writerow(["protein", "rhea_reaction"])
		for protein_id in mapping_prot_rhea:
			csvwriter.writerow([protein_id, ','.join(mapping_prot_rhea[protein_id])])
	print('\n')

	print('Find Uniprot protein with evidence for EC number')
	sparql.setQuery("""PREFIX up: <http://purl.uniprot.org/core/>
			PREFIX rh: <http://rdf.rhea-db.org/>
			PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

			SELECT DISTINCT ?protein ?enzyme WHERE {{
				# ECO 269 is experimental evidence
				BIND (<http://purl.obolibrary.org/obo/ECO_0000269> as ?evidence)

				?protein a up:Protein .
				?protein up:reviewed true.
				?protein up:enzyme ?enzyme .

				[] a rdf:Statement ;
						rdf:subject ?protein ;
						rdf:predicate up:enzyme ;
						rdf:object ?enzyme ;
						up:attribution ?attribution .
				?attribution up:evidence ?evidence .
			}}""")

	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	# Create mapping file uniprot ID -> Rhea Id with evidence.
	mapping_prot_ec = {}
	for result in results["results"]["bindings"]:
		protein_id = result["protein"]["value"].split('/')[-1]
		enzyme_id = result["enzyme"]["value"].split('/')[-1]
		if protein_id not in mapping_prot_ec:
			mapping_prot_ec[protein_id] = [enzyme_id]
		else:
			mapping_prot_ec[protein_id].append(enzyme_id)

	with open(database_folder + '/uniprot_ec_evidence.tsv', 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter ='\t')
		csvwriter.writerow(["protein", "EC_Number"])
		for protein_id in mapping_prot_ec:
			csvwriter.writerow([protein_id, ','.join(mapping_prot_ec[protein_id])])
	print('\n')

	print('Find Uniprot protein with evidence for GO term')
	sparql.setQuery("""PREFIX up: <http://purl.uniprot.org/core/>
	PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	SELECT DISTINCT ?protein ?goTerm
	WHERE
	{{
		# ECO 269 is experimental evidence
		BIND (<http://purl.obolibrary.org/obo/ECO_0000269> as ?evidence)

		?protein a up:Protein .
		?protein up:classifiedWith ?goTerm .
		Filter (contains(str(?goTerm),'http://purl.obolibrary.org/ob'))

		[] a rdf:Statement ;
				rdf:subject ?protein ;
				rdf:predicate up:classifiedWith ;
				rdf:object ?goTerm ;
				up:attribution ?attribution .
		?attribution up:evidence ?evidence .
	}}""")
	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	# Create mapping file uniprot ID -> Rhea Id with evidence.
	mapping_prot_go = {}
	for result in results["results"]["bindings"]:
		protein_id = result["protein"]["value"].split('/')[-1]
		goTerm_id = result["goTerm"]["value"].split('/')[-1]
		if protein_id not in mapping_prot_go:
			mapping_prot_go[protein_id] = [goTerm_id]
		else:
			mapping_prot_go[protein_id].append(goTerm_id)

	with open(database_folder + '/uniprot_go_evidence.tsv', 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter ='\t')
		csvwriter.writerow(["protein", "goTerm"])
		for protein_id in mapping_prot_go:
			csvwriter.writerow([protein_id, ','.join(mapping_prot_go[protein_id])])
	print('\n')

	annotated_proteins = []
	annotated_proteins.extend([prot for prot in mapping_prot_rhea])
	annotated_proteins.extend([prot for prot in mapping_prot_ec])
	annotated_proteins.extend([prot for prot in mapping_prot_go])
	annotated_proteins = set(annotated_proteins)
	print('Create fasta containing proteins with experimental evidence for Rhea reactions')
	records = []
	for record in SeqIO.parse(database_folder + '/uniprot_sprot.fasta', 'fasta'):
		record.id = record.id.split('|')[1]
		record.description = ''
		if record.id in annotated_proteins:
			records.append(record)
	SeqIO.write(records, database_folder+'/uniprot_rhea_evidence.fasta', 'fasta')
	print('\n')

	print('Create Rhea SBMl file')
	rhea_rdf_file = database_folder + '/rhea.rdf'
	uniprot_rhea_evidence = database_folder + '/uniprot_rhea_evidence.tsv'

	output_rhea_sbml = database_folder + '/rhea.sbml'
	output_rhea_variable_stoichiometric_reactions = database_folder + '/rhea_variable_stoich_reactions.tsv'
	rhea_to_sbml(rhea_rdf_file, uniprot_rhea_evidence, output_rhea_sbml, output_rhea_variable_stoichiometric_reactions)

	# Find metadata
	now = datetime.now()
	download_date = now.strftime('%d/%m/%Y')

	rhea_mapper_dependencies = ['cobra=='+cobra_version, 'biopython=='+biopython_version, 'rdflib=='+rdflib_version, 'sparqlwrapper=='+sparqlwrapper_version]

	# Create version file.
	rhea_mapper_version = pkg_resources.get_distribution("rhea_mapper").version
	versions = {'Download_date': download_date,
				'Rhea_release_number': rhea_release_number,
				'Rhea_release_date': rhea_release_date,
				'Swissprot_release_number': swissprot_release_number,
				'Swissprot_release_date': swissprot_release_date,
				'Trembl_release_number': trembl_release_number,
				'Trembl_release_date': trembl_release_date,
				'rhea_mapper': rhea_mapper_version,
				'rhea_mapper_dependencies': rhea_mapper_dependencies}
	with open(database_folder+'/version.json', 'w') as output_file:
		json.dump(versions, output_file, indent=4)
