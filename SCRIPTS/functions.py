__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

import pandas as pd
from os import path
import re
from snakemake.utils import validate
from snakemake.utils import min_version


configfile: "config.yaml"

MASTER_LIST = pd.read_table(config["MASTER_LIST"], dtype=str).set_index(["sample","lane"], drop=False)
MASTER_LIST.index = MASTER_LIST.index.set_levels([i.astype(str) for i in MASTER_LIST.index.levels])

# samples list from master_list

SAMPLES = [ sample for sample in MASTER_LIST["sample"]]
LANES = [ runid for runid in MASTER_LIST["lane"]]

# Global wildcard constraints

wildcard_constraints:
	sample = "|".join(MASTER_LIST["sample"]),
	lane = "|".join(MASTER_LIST["lane"]),



# # functions for start inputs rule

def fastqc_input(wildcards):
	""" if trimming use following dictionary"""
	if config["CUTADAPT"]["trimming"]:
		return {"R1": "FASTQ/TRIMMED/{sample}_{lane}_R1.fastq.gz".format(**wildcards),
		"R2":"FASTQ/TRIMMED/{sample}_{lane}_R2.fastq.gz".format(**wildcards)}
	else:
		"""
		function provides the rowwise fastqs from the master list in the form a python dictionary
		if no trimming use raw fastq files
		"""
		fastqs = MASTER_LIST.loc[(wildcards.sample, wildcards.lane), ["fastq1", "fastq2"]].dropna()
		if len(fastqs) == 2:
			return {"R1": fastqs.fastq1, "R2": fastqs.fastq2}
		else:
			return {"R1": fastqs.fastq1}

def cutadapt_input(wildcards):
		""" function provides the rowwise fastqs from the master list in the form a python dictionary"""
		fastqs = MASTER_LIST.loc[(wildcards.sample, wildcards.lane), ["fastq1", "fastq2"]].dropna()
		if len(fastqs) == 2:
			return {"R1": fastqs.fastq1, "R2": fastqs.fastq2}
		else:
			return {"R1": fastqs.fastq1}

def base_name(filepath):
	"""
	python function will get the files name after subtituting .fastq.gz extention from the fastq1 and fastq2 column of master list
	the file name will be used in the replacing the .fastq.gz with _fastq.html and _fastq.zip
	"""
	base = path.basename(filepath)
	base = re.sub("\\.fastq.gz$","",base)
	return base

def multiqc_input(wildcards):
	return expand(["QC/TRIMMED/{sample}_{lane}.paired.qc.txt","QC/FASTQC/{sample}_{lane}_R1_fastqc.zip"], zip,sample=SAMPLES, lane=LANES)

def salmon_input(wildcards):
	fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
	return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}


def star_input(wildcards):
	fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
	return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}