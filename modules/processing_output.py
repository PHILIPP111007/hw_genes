import os
import json
import multiprocessing
from datetime import datetime
from functools import partial

from .my_conf import (
	WORK_DIR,
	vcf_dir,
	samples_file,
	chrs_info_path,
	reference_dir,
	consensuses_dir,
	consensuses_file
)


def _get_chrs() -> list[str]:
	with open(chrs_info_path, 'r', encoding='UTF-8') as f:
		chrs: list[str] = json.load(f)

	return chrs


def _get_samples() -> list[str]:
	with open(samples_file, 'r', encoding='UTF-8') as f:
		samples: list[str] = json.load(f)

	return samples


def _get_vcf_data(chr_num: str) -> dict:
	vcf_file = f'{vcf_dir}/{chr_num}.json'
	with open(vcf_file, 'r', encoding='UTF-8') as f:
		vcf_data: dict = json.load(f)

	return vcf_data


def _get_reference_data(chr_num: str) -> str:
	reference_file = f'{reference_dir}/{chr_num}.txt'
	with open(reference_file, 'r', encoding='UTF-8') as f:
		dna = f.readline()

	return dna


def _get_consensuses() -> dict[str, list[int]]:
	with open(consensuses_file, 'r', encoding='UTF-8') as f:
		consensuses = json.load(f)

	return consensuses


def _create_fasta_for_sample(
		sample: str,
		consensuses_length: int,
		output_path: str,
		consensuses_per_chr
	) -> None:
	
	result_fasta_file = f'{output_path}/{sample}.fasta'
	with open(result_fasta_file, 'w', encoding='UTF-8') as f:

		for chr_num, val in consensuses_per_chr.items():
			for record in val:
				line = f'>chr{chr_num}:{record[0]}-{record[0] + consensuses_length}\n{record[1]}\n'
				f.write(line)


def _process_chr_num(
		chr_num: str,
		consensuses_per_chr,
		consensus_indexes_per_chr: dict[str, list[int]],
		consensuses_length: int,
		sample: str
	) -> None:
	
	print(f'chr{chr_num}')

	pos_min = min(consensus_indexes_per_chr[chr_num])
	pos_max = max(consensus_indexes_per_chr[chr_num]) + consensuses_length

	# Make DNA smaller with saved indexes for future calculations
	dna: list[str] = (_get_reference_data(chr_num=chr_num)[pos_min:pos_max]).split()
	vcf_data = _get_vcf_data(chr_num=chr_num)

	# Create set where DNA can change
	open_for_change = set()
	j = consensuses_length + 1 # промежуточное значение
	for n in consensus_indexes_per_chr[chr_num]:
		for i in range(n, n + j):
			open_for_change.add(i)

	for pos in vcf_data.keys():
		pos = int(pos)
		if pos not in open_for_change:
			continue

		dna = dna[:pos - pos_min] + [vcf_data[str(pos)]["alt"][sample]] + \
			dna[pos - pos_min + 1:]

	dna = ''.join(dna)
	consensuses_per_chr[chr_num] = [
		(
			i,
			dna[i - pos_min:i - pos_min + consensuses_length]
		) 
		for i in consensus_indexes_per_chr[chr_num]
	]


def prepare_dna_and_create_fasta_per_sample(
		consensuses_length: int,
		output_path: str
	) -> None:

	# Create folder for fasta files
	result_folder = f'/fasta_samples_{datetime.now()}'
	output_path += result_folder
	if not os.path.exists(output_path):
		os.mkdir(output_path)
	
	consensus_indexes_per_chr = _get_consensuses()
	chrs = _get_chrs()
	samples = _get_samples()

	for sample in samples:
		print(f'{sample = }: Start.')

		with multiprocessing.Manager() as manager:
			consensuses_per_chr = manager.dict()

			with multiprocessing.Pool() as pool:

				partial_process_chr_num = partial(
					_process_chr_num, 
					consensuses_per_chr=consensuses_per_chr, 
					consensus_indexes_per_chr=consensus_indexes_per_chr, 
					consensuses_length=consensuses_length, 
					sample=sample
				)

				pool.map_async(partial_process_chr_num, chrs)
				pool.close()
				pool.join()

			_create_fasta_for_sample(
				sample=sample,
				consensuses_length=consensuses_length,
				output_path=output_path,
				consensuses_per_chr=consensuses_per_chr
			)
			print(f'{sample = }: Done.')


def delete_work_folders() -> None:
	"""Deletes work folders."""

	for directory_path in (vcf_dir, reference_dir, consensuses_dir, WORK_DIR):
		try:
			paths = os.listdir(directory_path)
			for path in paths:
				path = os.path.join(directory_path, path)
				if os.path.isfile(path):
					os.remove(path)
				else:
					os.rmdir(path)
		except Exception as e:
			print(e)

	os.rmdir(WORK_DIR)
