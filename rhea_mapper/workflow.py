import os

from rhea_mapper import create_database, input_parser, rhea_reconstruction

def rhea_mapper_workflow(input_folder, output_folder, database_folder):
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
		rhea_reconstruction.manage_genome(input_folder, input_tsv, input_fasta, database_folder, output_folder)