import os
import csv
import json
import multiprocessing
from itertools import groupby

from .my_conf import (
	WORK_DIR,
	vcf_dir,
	samples_file,
	chrs_info_path,
	reference_dir,
	consensuses_dir
)


def _prepare_vcf(vcf_path: str) -> str:
	file_new_path = vcf_dir + '/' + vcf_path.split('/')[-1]
	
	with open(vcf_path, 'r', encoding='UTF-8') as old_file:
		data = old_file.readlines()

	i = 0
	for _ in range(len(data)):
		if data[i].startswith('##'):
			data.pop(0)
		else:
			i += 1

	with open(file_new_path, 'w', encoding='UTF-8') as new_file:
		new_file.writelines(data)

	return file_new_path


def _get_samples(vcf_path: str) -> list[str]:
	"""Get list of samples and save those to a new file."""

	with open(vcf_path, 'r', encoding='UTF-8') as f:
		reader = csv.DictReader(f, delimiter='\t')
		samples: list[str] = list(reader.__next__().keys())[9:]

	with open(samples_file, 'w', encoding='UTF-8') as f:
		json.dump(samples, f)

	return samples


def _get_fasta_iter(fasta_path: str) -> tuple[str, str]:
	"""
	Given a fasta file. yield tuples of header, sequence.
	"""

	with open(fasta_path, 'r', encoding='UTF-8') as f:	
		# ditch the boolean (x[0]) and just keep the header or sequence since
		# we know they alternate.
		faiter = (x[1] for x in groupby(f, lambda line: line[0] == '>'))
		for header in faiter:
			# drop the '>'
			header = header.__next__().strip()
			if header.startswith('>chr'):
				header = header[4:]
			elif header.startswith('>'):
				header = header[1:]

			# join all sequence lines to one.
			seq = ''.join(s.strip() for s in faiter.__next__())
			yield header, seq


def _create_fasta_file(line, chrs):
	header = line[0]
	if header in chrs:
		print(f'chr{header}')

		file_path = f'{reference_dir}/{header}.txt'
		with open(file_path, 'w', encoding='UTF-8') as f:
			f.write(line[1])


def create_work_folders() -> None:
	"""Creates work folders."""

	for path in (WORK_DIR, reference_dir, vcf_dir, consensuses_dir):
		if not os.path.exists(path):
			os.mkdir(path)


def vcf_processing(vcf_path: str, AF: float) -> None:
	"""
	Reads VCF file with multiple samples.
	Save in the new file a list with chromosomes names.
	Save prepared info about each chromosome.
	"""

	vcf_path = _prepare_vcf(vcf_path=vcf_path)
	samples = _get_samples(vcf_path=vcf_path)
	DB_vcf = {}
	chrs = set()

	with open(vcf_path, 'r', encoding='UTF-8') as f:
		reader = csv.DictReader(f, delimiter='\t')

		for record in reader:
			chr_num = record['#CHROM']

			if chr_num.startswith('chr'):
				chr_num = chr_num[3:]

			pos = record['POS']
			chrs.add(chr_num)

			# Include record based on AF.
			include_record = True
			for kv in record['INFO'].split(';'):
				kv = kv.split('=')
				if kv[0] == 'AF':
					af_lst = kv[1].split(',')
					for af in af_lst:
						if float(af) < AF:
							include_record = False
							break
					break

			if not include_record:
				continue

			result = {} # {'SAMPLE1': 'GT', 'SAMPLE2': 'AT'}
			alleles: list[str] = list(record['REF']) + record['ALT'].split(',')
			delimeter = '/' # SAMPLE1 - 0/1
			for sample in samples:
				allele_indexes = record[sample].split(':')[0]
				if '|' in allele_indexes:
					delimeter = '|'

				allele_indexes = allele_indexes.split(delimeter)

				if ''.join(allele_indexes).isdigit():
					allele_indexes = list(map(int, allele_indexes))
					result[sample] = ''.join(list(alleles[i] for i in allele_indexes))

			if result:
				DB_vcf.setdefault(chr_num, {})
				DB_vcf[chr_num].setdefault(pos, result)

	chrs: list[str] = list(chrs)

	# Save a list with chromosomes names
	with open(chrs_info_path, 'w', encoding='UTF-8') as f:
		json.dump(list(chrs), f)

	# Save prepared info about each chromosome
	for chr_num, val in DB_vcf.items():
		file_path = f'{vcf_dir}/{chr_num}.json'

		with open(file_path, 'w', encoding='UTF-8') as vcf_new_file:
			json.dump(val, vcf_new_file)


def fasta_processing(fasta_path: str) -> None:
	"""
	Creates new files, each file contains only one chromosome sequencing.
	"""

	fasta_iter = _get_fasta_iter(fasta_path=fasta_path)

	with open(chrs_info_path, 'r', encoding='UTF-8') as f:
		chrs: list[str] = json.load(f)

	with multiprocessing.Pool() as pool:
		pool.starmap(_create_fasta_file, [(line, chrs) for line in fasta_iter])
		pool.close()
		pool.join()
