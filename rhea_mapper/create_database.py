import csv
import gzip
import os
import sys
import shutil
import urllib.request

from Bio import SeqIO
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model
from rdflib import Graph
from SPARQLWrapper import SPARQLWrapper, JSON


def urllib_reporthook(count, block_size, total_size):
	downloaded_size = int(count * block_size) / (1024 * 1024)
	percent = min(int(count*block_size*100/total_size),100)
	total_size_show = total_size / (1024 * 1024)
	sys.stdout.write("\rDownloading: {0}%, {1:.2f} MB on {2:.2f} MB".format(percent, downloaded_size, total_size_show))
	sys.stdout.flush()


def rhea_to_sbml(rhea_rdf_file, uniprot_rhea_evidence, output_file):
	# Read rhea rdf file.
	g = Graph()
	g.parse(rhea_rdf_file)

	# Get all the approved reactions with their metabolites.
	query_reactions_metabolites = g.query(
		"""PREFIX rh:<http://rdf.rhea-db.org/>
		PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
		SELECT ?reaction ?reactionSide ?compound ?compoundName ?compoundFormula ?coefficient WHERE {
		?reaction rdfs:subClassOf rh:Reaction .
		?reaction rh:status rh:Approved .
		?reaction rh:side ?reactionSide .
		?reactionSide rh:contains ?participant .
		?reactionSide ?rhContains ?participant .
		?rhContains rdfs:subPropertyOf rh:contains .
		?rhContains rh:coefficient ?coefficient .
		?participant rh:compound ?compound .
		?compound rh:name ?compoundName .
		?compound rh:formula ?compoundFormula .
		}""")

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
		SELECT ?reaction ?bidirection ?xref WHERE {
		?reaction rdfs:subClassOf rh:Reaction .
		?reaction rh:status rh:Approved .
		?reaction rh:side ?reactionSide .
		?reaction rh:bidirectionalReaction ?bidirection .
		?bidirection rdfs:seeAlso ?xref .
		}""")

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
		SELECT ?reaction ?compoundproduct ?compoundsubstrate ?xref WHERE {
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
		}""")

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
	for row in query_reactions_metabolites:
		reaction_id, reaction_side, metabolite_id, metabolite_name, metabolite_formula, stoichiometric_coefficient = row
		if reaction_id not in reactions:
			reaction_id_sbml = reaction_id.split('/')[-1]
			metabolite_id = metabolite_id.split('/')[-1]
			# Create metabolite.
			if metabolite_id not in metabolites:
				metabolite_id_sbml = Metabolite(metabolite_id, compartment='c', name=metabolite_name, formula=metabolite_formula)
				metabolites[metabolite_id] = metabolite_id_sbml
			else:
				metabolite_id_sbml = metabolites[metabolite_id]
			if reaction_id_sbml not in reactions:
				reactions[reaction_id_sbml] = {}

			# Manage specific stoechiometric coefficient.
			if str(stoichiometric_coefficient) == "N":
				stoichiometric_coefficient = 2.0
			elif str(stoichiometric_coefficient) == "Nplus1":
				stoichiometric_coefficient = 3.0
			elif str(stoichiometric_coefficient) == "Nminus1":
				stoichiometric_coefficient = 1.0
			elif str(stoichiometric_coefficient) == "2n":
				stoichiometric_coefficient = 4.0
			else:
				stoichiometric_coefficient = float(stoichiometric_coefficient)

			# Create reation in the direction find by the directed reaction SPARQL query.
			if reaction_id_sbml in directed_reactions:
				# Product compounds are assigned to 1.0 in cobra model.
				if metabolite_id in directed_reactions[reaction_id_sbml]["product"]:
					reactions[reaction_id_sbml][metabolite_id_sbml] = stoichiometric_coefficient
				# Substrate compounds are assigned to -1.0 in cobra model.
				elif metabolite_id in directed_reactions[reaction_id_sbml]["substrate"]:
					reactions[reaction_id_sbml][metabolite_id_sbml] = - stoichiometric_coefficient
			#Otherwise use the arbitrary direction.
			else:
				# Right compounds are assigned to 1.0 in cobra model.
				if reaction_side.endswith('_R'):
					reactions[reaction_id_sbml][metabolite_id_sbml] = stoichiometric_coefficient
				# Left compounds are assigned to -1.0 in cobra model.
				elif reaction_side.endswith('_L'):
					reactions[reaction_id_sbml][metabolite_id_sbml] = - stoichiometric_coefficient

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
	write_sbml_model(model, output_file)


def download_database(database_folder):
	if not os.path.exists(database_folder):
		os.mkdir(database_folder)

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

	print('Find Uniprot protein with experimental evidence for Rhea reaction')
	sparql = SPARQLWrapper('https://sparql.uniprot.org/sparql')

	sparql.setQuery("""PREFIX up: <http://purl.uniprot.org/core/>
			PREFIX rh: <http://rdf.rhea-db.org/>
			PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

			SELECT distinct ?rhea_reaction ?protein ?protein_name WHERE {

			# ECO 269 is experimental evidence
			BIND (<http://purl.obolibrary.org/obo/ECO_0000269> as ?evidence)
			?protein up:recommendedName/up:fullName ?protein_name .
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
			}""")

	# Parse output.
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()

	# Create mapping file uniprot ID -> Rhea Id with evidence.
	mapping_prot_rhea = {}
	for result in results["results"]["bindings"]:
		protein_id = result["protein"]["value"].split('/')[-1]
		reaction_id = result["rhea_reaction"]["value"].split('/')[-1]
		protein_name = result["protein_name"]["value"].split('/')[-1]
		if protein_id not in mapping_prot_rhea:
			mapping_prot_rhea[protein_id] = ([reaction_id], protein_name)
		else:
			if reaction_id not in mapping_prot_rhea[protein_id][0]:
				mapping_prot_rhea[protein_id][0].append(reaction_id)

	with open(database_folder + '/uniprot_rhea_evidence.tsv', 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter ='\t')
		csvwriter.writerow(["protein", "rhea_reaction", "protein_name"])
		for protein_id in mapping_prot_rhea:
			reaction_ids = mapping_prot_rhea[protein_id][0]
			protein_name = mapping_prot_rhea[protein_id][1]
			csvwriter.writerow([protein_id, ','.join(reaction_ids), protein_name])
	print('\n')

	print('Create fasta containing proteins with experimental evidence for Rhea reactions')
	records = []
	for record in SeqIO.parse(database_folder + '/uniprot_sprot.fasta', 'fasta'):
		record.id = record.id.split('|')[1]
		record.description = ''
		if record.id in mapping_prot_rhea:
			records.append(record)
	SeqIO.write(records, database_folder+'/uniprot_rhea_evidence.fasta', 'fasta')
	print('\n')

	print('Create Rhea SBMl file')
	rhea_to_sbml(database_folder + '/rhea.rdf', database_folder + '/uniprot_rhea_evidence.tsv', database_folder + '/rhea.sbml')
