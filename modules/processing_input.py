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

	for _ in range(len(data)):
		if data[0].startswith('##'):
			data.pop(0)

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
			else:
				header = header[1:]

			# join all sequence lines to one.
			seq = ''.join(s.strip().upper() for s in faiter.__next__())
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

	if not os.path.exists(WORK_DIR):
		os.mkdir(WORK_DIR)

	if not os.path.exists(reference_dir):
		os.mkdir(reference_dir)

	if not os.path.exists(vcf_dir):
		os.mkdir(vcf_dir)

	if not os.path.exists(consensuses_dir):
		os.mkdir(consensuses_dir)


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
			alt = record['ALT'].split(',')

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

			result = {
				samples[i]: alt[i] if len(alt) > 1 else alt[0]
				for i in range(len(samples))
			} # {'SAMPLE1': 'G', 'SAMPLE2': 'A'}

			DB_vcf.setdefault(chr_num, {})
			DB_vcf[chr_num].setdefault(pos, {
				'ref': record['REF'],
				'alt': {}
			})['alt'] = result

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
