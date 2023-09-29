#!/usr/local/bin/python3

import argparse
import time


from modules.processing_input import (
	create_work_folders,
	fasta_processing,
	vcf_processing
)
from modules.processing_middle import (
	create_consensus_indexes
)
from modules.processing_output import (
	prepare_dna_and_create_fasta_per_sample,
	delete_work_folders
)


parser = argparse.ArgumentParser()
parser.add_argument('-r', '--reference', required=True, help='Reference file.')
parser.add_argument('-vcf', '--vcf', required=True, help='VCF file.')
parser.add_argument('-o', '--output', required=True, help='Output directory.')
parser.add_argument('-l', '--length', type=int, default=100, 
					help='Consensuses length (default is 100).')
parser.add_argument('-c', '--count', type=int, default=1_000_000, 
					help='Consensuses count (default is 1_000_000).')
parser.add_argument('-af', '--allele_frequency', type=float, default=0.5, 
					help='Allele frequency (default is 0.5).')

args = parser.parse_args()


if __name__ == '__main__':

	fasta_path: str = args.reference
	vcf_path: str = args.vcf
	output_path: str = args.output
	consensuses_length: int = args.length
	consensuses_count: int = args.count
	AF: float = args.allele_frequency

	t_1 = time.time()

	create_work_folders()

	print('[1 / 4] VCF file: Start.')
	vcf_processing(vcf_path=vcf_path, AF=AF)
	print('[1 / 4] VCF file: Done.')

	print('[2 / 4] Fasta file: Start.')
	fasta_processing(fasta_path=fasta_path)
	print('[2 / 4] Fasta file: Done.')

	print('[3 / 4] Consensus indexes: Start.')
	create_consensus_indexes(
		consensuses_length=consensuses_length,
		consensuses_count=consensuses_count
	)
	print('[3 / 4] Consensus indexes: Done.')

	print('[4 / 4] prepare_dna_and_create_fasta_per_sample: Start.')
	prepare_dna_and_create_fasta_per_sample(
		consensuses_length=consensuses_length,
		output_path=output_path
	)
	print('[4 / 4] prepare_dna_and_create_fasta_per_sample: Done.')

	delete_work_folders()

	t_2 = time.time()
	z = (t_2 - t_1) / 60
	print(f'Time (min): {z}')
