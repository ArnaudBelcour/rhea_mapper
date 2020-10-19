import csv
import subprocess
import shutil
import os

from rdflib import Graph
from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model, write_sbml_model

def read_tsv(tsv_pathname, splitter=None):
	tsv_dict = {}
	with open(tsv_pathname, 'r') as tsv_file:
		csvreader = csv.reader(tsv_file, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			if splitter:
				tsv_dict[line[0]] = line[1].split(splitter)
			else:
				tsv_dict[line[0]] = line[1]
	
	return tsv_dict

def map_ec_to_rhea(annotation_pathname, rhea_ec_pathname, output_pathname):
	gene_ecs = read_tsv(annotation_pathname, splitter=',')

	rhea_ecs = {}
	with open(rhea_ec_pathname, 'r') as tsv_file:
		csvreader = csv.reader(tsv_file, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			if len(line) > 0:
				rhea_ecs[line[3]] = line[0]


	gene_rheas = {}
	for gene_id in gene_ecs:
		gene_rheas[gene_id] = [rhea_ecs[ec] for ec in gene_ecs[gene_id] if ec in rhea_ecs]

	with open(output_pathname, 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter='\t')
		csvwriter.writerow(['gene', 'Rhea'])
		for gene in gene_rheas:
			csvwriter.writerow([gene, ','.join(gene_rheas[gene])])

def orthology_to_uniprot(input_proteome, rhea_proteome, output_folder, orthogroups_result, nb_cpu):
	work_directory = output_folder + '/orthology/'
	if os.path.exists(work_directory):
		shutil.rmtree(work_directory)
	os.mkdir(work_directory)
	input_proteome_file = os.path.basename(input_proteome)
	shutil.copyfile(input_proteome, work_directory+input_proteome_file)
	rhea_proteome_file = os.path.basename(rhea_proteome)
	shutil.copyfile(rhea_proteome, work_directory+rhea_proteome_file)

	subprocess.call(['orthofinder', '-f', work_directory, '-og', '-t', str(nb_cpu)])

	orthodata_path = max(["%s/%s" %(x[0], 'Orthogroups/Orthogroups.tsv') for x in os.walk(work_directory) if 'Orthogroups' in x[1]])
	shutil.copyfile(orthodata_path, orthogroups_result)


def map_orthology_to_rhea(rhea_uniprot_pathname, orthofinder_result_pathname, output_pathname):
	uniprot_to_rheas = {}
	with open(rhea_uniprot_pathname, 'r') as rhea_mapping:
		csvreader = csv.reader(rhea_mapping, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			reaction_ids = line[1].split(',')
			protein_id = line[0]
			if protein_id not in uniprot_to_rheas:
				uniprot_to_rheas[protein_id] = reaction_ids
			else:
				uniprot_to_rheas[protein_id].extend(reaction_ids)

	gene_rheas = {}  

	with open(orthofinder_result_pathname, 'r') as orthofinde_results:
		csvreader = csv.reader(orthofinde_results, delimiter='\t')
		next(csvreader)
		for line in csvreader:
			if line[1] != '' and line[2] != '':
				gene_prots = [gene_id.split("_isoform")[0] for gene_id in line[1].split(", ")]
				uniprot_prots = [gene_id for gene_id in line[2].split(", ")]
				rhea_reactions = list(set([rhea_reaction
								for uniprot_prot in uniprot_prots if uniprot_prot in uniprot_to_rheas
								for rhea_reaction in uniprot_to_rheas[uniprot_prot]]))
				if gene_prots != [] and uniprot_prots != []:
					for gene_prot in gene_prots:
						gene_rheas[gene_prot] = rhea_reactions

	with open(output_pathname, 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter='\t')
		csvwriter.writerow(['gene', 'Rhea'])
		for gene in gene_rheas:
			csvwriter.writerow([gene, ','.join(gene_rheas[gene])])


def merge_mapping_orthology(mapping_pathname, orthology_pathname, output_file):
	results = {}

	mapping_rheas = read_tsv(mapping_pathname, splitter=',')

	orthology_rheas = read_tsv(orthology_pathname, splitter=',')

	all_genes = list(mapping_rheas.keys()) + list(orthology_rheas.keys())

	for gene in all_genes:
		if gene in mapping_rheas:
			gene_mapping_rheas = mapping_rheas[gene]
		else:
			gene_mapping_rheas = []
		if gene in orthology_rheas:
			gene_orthology_rheas = orthology_rheas[gene]
		else:
			gene_orthology_rheas = []
		results[gene] = list(set(gene_mapping_rheas).union(gene_orthology_rheas))

	with open(output_file, 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter='\t')
		csvwriter.writerow(['gene', 'Rhea'])
		for gene in results:
			csvwriter.writerow([gene, ','.join(results[gene])])


def sbml_creation(rhea_sbml_file, rhea_mapping_file, sbml_output_file):
	rhea_model = read_sbml_model(rhea_sbml_file)

	genome_reactions = {}
	with open(rhea_mapping_file, 'r') as rhea_mapping:
		csvreader = csv.reader(rhea_mapping, delimiter='\t')
		for line in csvreader:
			for reaction in line[1].split(','):
				if reaction not in genome_reactions:
					genome_reactions[reaction] = [line[0]]
				else:
					genome_reactions[reaction].append(line[0])

	species_model = Model('test')

	sbml_reactions = []
	for reaction in rhea_model.reactions:
		if reaction.id.replace('R_','') in genome_reactions:
			reaction.gene_reaction_rule = '( ' + ' or '.join([prot for prot in genome_reactions[reaction.id.replace('R_','')]]) + ' )'
			sbml_reactions.append(reaction)

	species_model.add_reactions(sbml_reactions)

	# Create sbml file.
	write_sbml_model(species_model, sbml_output_file)

def manage_genome(input_folder, annotation_pathname, input_proteome, database_folder, output_folder):
	rhea2ec = database_folder + '/rhea2ec.tsv'
	uniprot_rhea_evidence_fasta = database_folder + '/uniprot_rhea_evidence.fasta'
	uniprot_rhea_evidence_tsv = database_folder + '/uniprot_rhea_evidence.tsv'
	rhea_sbml = database_folder + '/rhea.sbml'

	output_tmp_folder = output_folder + '/tmp/' + input_folder

	if not os.path.exists(output_tmp_folder):
		os.mkdir(output_tmp_folder)

	mapping_result = output_tmp_folder + '/mapping.tsv'
	orthogroups_result = output_tmp_folder + '/orthogroups_result.tsv'
	orthology_result = output_tmp_folder +  '/orthology.tsv'
	merged_result = output_tmp_folder +  '/merged.tsv'
	metabolic_network = output_folder + '/' + input_folder + '.sbml'

	map_ec_to_rhea(annotation_pathname, rhea2ec, mapping_result)
	orthology_to_uniprot(input_proteome, uniprot_rhea_evidence_fasta, output_tmp_folder, orthogroups_result, 3)
	map_orthology_to_rhea(uniprot_rhea_evidence_tsv, orthogroups_result, orthology_result)
	merge_mapping_orthology(mapping_result, orthology_result, merged_result)
	sbml_creation(rhea_sbml, merged_result, metabolic_network)