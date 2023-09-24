import os
import json
import random
import multiprocessing

from .my_conf import (
	vcf_dir,
	chrs_info_path,
	reference_dir,
	consensuses_file
)


def _get_chrs() -> list[str]:
	with open(chrs_info_path, 'r', encoding='UTF-8') as f:
		chrs: list[str] = json.load(f)
	return chrs


def _check_files_exists(chr_num: str) -> bool:
	vcf_file = f'{vcf_dir}/{chr_num}.json'
	reference_file = f'{reference_dir}/{chr_num}.txt'

	return os.path.exists(vcf_file) and os.path.exists(reference_file)


def _get_vcf_data(chr_num: str) -> dict:
	vcf_file = f'{vcf_dir}/{chr_num}.json'
	with open(vcf_file, 'r', encoding='UTF-8') as f:
		vcf_data: dict = json.load(f)

	return vcf_data


def _get_distribution_across_chromosomes(
		consensuses_count: int,
		chrs: list[str]
	) -> dict[str, int]:
	
	number_per_chr = {}
	N = consensuses_count

	n = N // len(chrs)

	for chr_num in chrs:
		if not _check_files_exists(chr_num=chr_num):
			continue
		number_per_chr[chr_num] = n

	number_per_chr[chrs[0]] += N % len(chrs)
	return number_per_chr


def _get_random_distribution_across_one_chromosome(
		vcf_data: dict,
		count: int
	) -> dict[str, int]:

	positions: list[str] = list(vcf_data.keys())
	distribution = [random.choice(positions) for _ in range(count)]
	result = {pos: distribution.count(pos) for pos in positions}
	return result


def _save_consensuses(consensuses) -> None:
	with open(consensuses_file, 'w', encoding='UTF-8') as f:
		json.dump(consensuses, f)


def _process_chr_num(
		consensuses_per_chr,
		consensuses_length: int,
		chr_num: str,
		count: int
	):

	print(f'chr{chr_num}')
	vcf_data = _get_vcf_data(chr_num=chr_num)
	distribution: dict[str: int] = _get_random_distribution_across_one_chromosome(
		vcf_data=vcf_data,
		count=count
	)

	consensuses = []
	i = consensuses_length + 1 # промежуточное значение
	for pos, count in distribution.items():
		pos = int(pos)
		j = pos - i
		if j < 0:
			j = 0

		for _ in range(count):
			"""Get random start position on the DNA."""

			consensuses.append(random.randint(j, pos))
	consensuses = set(consensuses)
	consensuses_per_chr[chr_num] = sorted(consensuses)


def create_consensus_indexes(consensuses_length: int, consensuses_count: int):
	chrs: list = _get_chrs()
	chrs: dict[str, int] = _get_distribution_across_chromosomes(
		consensuses_count=consensuses_count,
		chrs=chrs
	)

	# Update chromosomes
	with open(chrs_info_path, 'w', encoding='UTF-8') as f:
		json.dump(list(chrs.keys()), f)

	with multiprocessing.Manager() as manager:
		consensuses_per_chr = manager.dict()

		with multiprocessing.Pool() as pool:

			pool.starmap(_process_chr_num, [(
				consensuses_per_chr,
				consensuses_length,
				chr_num,
				count
			) for chr_num, count in chrs.items()])

			pool.close()
			pool.join()

		_save_consensuses(consensuses=consensuses_per_chr.copy())
