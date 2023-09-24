# File contains constants.

from pathlib import Path


WORK_DIR: str = Path(__file__).resolve().parent.parent.as_posix() + '/work_dir'

vcf_dir: str = f'{WORK_DIR}/vcf'
chrs_info_path: str = f'{vcf_dir}/chrs.json'
samples_file: str = f'{vcf_dir}/samples.json'
reference_dir: str = f'{WORK_DIR}/reference'
consensuses_dir: str = f'{WORK_DIR}/consensuses'
consensuses_file: str = f'{consensuses_dir}/consensuses.json'
