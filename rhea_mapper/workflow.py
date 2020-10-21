import os

from rhea_mapper import create_database, input_parser, rhea_reconstruction, sbml_from_file, sparql_query


def rhea_mapper_workflow(input_folder, output_folder, database_folder, nb_cpu=1):
	if not os.path.exists(database_folder):
		print('Missing database folder, it will be created in ' + database_folder)
		create_database.download_database(database_folder)

	if not os.path.exists(output_folder):
		print('Create output folder at ' + output_folder)
		os.mkdir(output_folder)
		os.mkdir(output_folder + '/tmp')

	rhea_mapper_input_data = input_parser.check_input(input_folder)

	for input_folder in rhea_mapper_input_data:
		input_fasta, input_tsv = rhea_mapper_input_data[input_folder]
		rhea_reconstruction.manage_genome(input_folder, input_tsv, input_fasta, database_folder, output_folder, nb_cpu)


def sparql_query_workflow(user_sparql_query, output_folder, database_folder, endpoint, nb_cpu):
	if not os.path.exists(output_folder):
		print('Create output folder at ' + output_folder)
		os.mkdir(output_folder)

	if endpoint == 'uniprot':
		sparql_query.query_uniprot(user_sparql_query, database_folder, output_folder, nb_cpu)
	elif endpoint == 'rhea':
		sparql_query.query_rhea(user_sparql_query, database_folder, output_folder)
		sparql_query.rhea_sbml_creation(database_folder+'/rhea.sbml', output_folder)
	elif endpoint == 'rhea_endpoint':
		sparql_query.query_rhea_endpoint(user_sparql_query, output_folder)
		sparql_query.rhea_sbml_creation(database_folder+'/rhea.sbml', output_folder)


def sbml_from_list(element_pathname, database_folder, output_folder, database, nb_cpu):
	if not os.path.exists(output_folder):
		print('Create output folder at ' + output_folder)
		os.mkdir(output_folder)

	if database == 'uniprot':
		sbml_from_file.uniprot_from_file(element_pathname, database_folder+'/rhea.sbml', database_folder, output_folder, nb_cpu)
	elif database == 'rhea':
		sbml_from_file.rhea_sbml_from_file(element_pathname, database_folder+'/rhea.sbml', output_folder)
