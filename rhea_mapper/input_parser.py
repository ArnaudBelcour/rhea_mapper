import csv
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    # Import to be compatible with biopython version lesser than 1.78
    from Bio.Alphabet.IUPAC import protein
except ImportError:
    # Exception to be compatible with biopython version superior to 1.78
    protein = None


def genbank_parser(genbank_file, fasta_output, annotation_output):
	fasta_records = []
	fasta_ids = {}
	fastas_ecs = {}
	with open(genbank_file, 'r') as gbk:
		for seq_record in SeqIO.parse(gbk, 'genbank'):
			seq_feature_cds = (seq_feature for seq_feature in seq_record.features if seq_feature.type == "CDS")
			for seq_feature in seq_feature_cds:
				if 'locus_tag' in seq_feature.qualifiers:
					fasta_id = seq_feature.qualifiers['locus_tag'][0]
					if 'EC_number' in seq_feature.qualifiers:
						if fasta_id not in fastas_ecs:
							fastas_ecs[fasta_id] = seq_feature.qualifiers['EC_number']
						else:
							fastas_ecs[fasta_id].append(seq_feature.qualifiers['EC_number'])
					if 'translation' in seq_feature.qualifiers:
						if fasta_id not in fasta_ids:
							# Keep compatibility with biopython version lesser than 1.78
							if protein:
								fasta_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0], protein), id=fasta_id, description=fasta_id)
							else:
								fasta_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0]), id=fasta_id, description=fasta_id)
							fasta_records.append(fasta_record)
							fasta_ids[fasta_id] = 1
						else:
							fasta_ids[fasta_id] += 1
							isoform_id = fasta_id + '_isoform' + str(fasta_ids[fasta_id])	
							# Keep compatibility with biopython version lesser than 1.78
							if protein:
								fasta_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0], protein), id=isoform_id, description=fasta_id)
							else:
								fasta_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0]), id=isoform_id, description=fasta_id)
							fasta_records.append(fasta_record)

	SeqIO.write(fasta_records, fasta_output, "fasta")

	with open(annotation_output, 'w') as output_file:
		csvwriter = csv.writer(output_file, delimiter='\t')
		csvwriter.writerow(['gene', 'EC_Number'])
		for gene in fastas_ecs:
			csvwriter.writerow([gene, ','.join(fastas_ecs[gene])])

def check_input(input_folder_pathname):

	rhea_mapper_input_data = {}
	for input_folder in os.listdir(input_folder_pathname):
		input_genbank = input_folder_pathname + '/' + input_folder + '/' + input_folder + '.gbk'
		input_fasta = input_folder_pathname + '/' + input_folder + '/' + input_folder + '.fasta'
		input_tsv = input_folder_pathname + '/' + input_folder + '/' + input_folder + '.tsv'
		if os.path.exists(input_genbank):
			genbank_parser(input_genbank, input_fasta, input_tsv)

		if not os.path.exists(input_fasta):
			print('Missing fasta file for ' + input_folder)
		if not os.path.exists(input_tsv):
			print('Missing tsv file for ' + input_folder)
		rhea_mapper_input_data[input_folder] = (input_fasta, input_tsv)

	return rhea_mapper_input_data