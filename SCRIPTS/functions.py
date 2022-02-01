__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

configfile: "config.yaml"

MASTER_LIST = pd.read_table(config["MASTER_LIST"], dtype=str).set_index(["sample"], drop=False)
#MASTER_LIST.index = MASTER_LIST.index.set_levels([i.astype(str) for i in MASTER_LIST.index.levels])

# samples list from master_list

SAMPLES = [ sample for sample in MASTER_LIST["sample"]]

# wildcard constraints

wildcard_constraints:
	sample = "|".join(MASTER_LIST["sample"])



# # functions for start inputs rule

def star_input(wildcards):
	fastqs=MASTER_LIST.loc[(wildcards.sample),["fastq1","fastq2"]].dropna()
	for f in [fastqs.fastq1]:
		if len(fastqs.fastq1) < 5:
			return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}
		else:
			return {"R1": fastqs.fastq1, "R2": fastqs.fastq1}


def salmon_input(wildcards):
	fastqs=MASTER_LIST.loc[(wildcards.sample),["fastq1","fastq2"]].dropna()
	for f in [fastqs.fastq1]:
		if len(fastqs.fastq1) < 5:
			return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}
		else:
			return {"R1": fastqs.fastq1, "R2": fastqs.fastq1}



# def star_input(wildcards):
# 	""" function provides the rowwise fastqs from the master list in the form a python dictionary"""
# 	fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
# 	if len(fastqs)> 2:
# 		return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}
# 	elif len(fastqs)==2:
# 		return {"R1": fastqs.fastq1, "R2": fastqs.fastq2}
# 	else:
# 		return {"R1": [f for f in fastqs.fastq1]}

# def salmon_input(wildcards):
# 	""" function provides the rowwise fastqs from the master list in the form a python dictionary"""
# 	fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
# 	if len(fastqs) > 2:
# 		return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}
# 	elif len(fastqs)==2:
# 		return {"R1": fastqs.fastq1, "R2": fastqs.fastq2}
# 	else:
# 		return {"R1": fastqs.fastq1}