__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_DNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

import pandas as pd
from os import path
import re
from snakemake.utils import validate
from snakemake.utils import min_version


configfile: "config.yaml"

MASTER_LIST = pd.read_table(config["MASTER_LIST"], dtype=str).set_index(["sample","runID"], drop=False)
MASTER_LIST.index = MASTER_LIST.index.set_levels([i.astype(str) for i in MASTER_LIST.index.levels])

# samples list from master_list

SAMPLES = [ sample for sample in MASTER_LIST["sample"]]
RUNIDS = [ runID for runID in MASTER_LIST["runID"]]

# Global wildcard constraints

wildcard_constraints:
	sample = "|".join(MASTER_LIST["sample"]),
	runID = "|".join(MASTER_LIST["runID"]),



# # functions for start inputs rule

def fastqc_input(wildcards):
	""" if trimming use following dictionary"""
	if config["CUTADAPT"]["trimming"]:
		return {"R1": "FASTQ/TRIMMED/{sample}_{runID}_R1.fastq.gz".format(**wildcards),
		"R2":"FASTQ/TRIMMED/{sample}_{runID}_R2.fastq.gz".format(**wildcards)}
	else:
		"""
		function provides the rowwise fastqs from the master list in the form a python dictionary
		if no trimming use raw fastq files
		"""
		fastqs = MASTER_LIST.loc[(wildcards.sample, wildcards.runID), ["fastq1", "fastq2"]].dropna()
		if len(fastqs) == 2:
			return {"R1": fastqs.fastq1, "R2": fastqs.fastq2}
		else:
			return {"R1": fastqs.fastq1}

def cutadapt_input(wildcards):
		""" function provides the rowwise fastqs from the master list in the form a python dictionary"""
		fastqs = MASTER_LIST.loc[(wildcards.sample, wildcards.runID), ["fastq1", "fastq2"]].dropna()
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
	if config["CUTADAPT"]["trimming"]:
		return expand(["QC/TRIMMED/{sample}_{runID}.paired.qc.txt",
			"QC/FASTQC/{sample}_{runID}_R1_fastqc.zip",
			"QC/FASTQC/{sample}_{runID}_R2_fastqc.zip"],zip, sample=SAMPLES, runID=RUNIDS)
	else:
		return expand(["QC/FASTQC/{sample}_{runID}_R1_fastqc.zip",
			"QC/FASTQC/{sample}_{runID}_R2_fastqc.zip"],zip, sample=SAMPLES, runID=RUNIDS)

def salmon_input(wildcards):
	if not config["CUTADAPT"]["trimming"]:
		fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
		return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}
	else:
		fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
		runs = fastqs.index.tolist()
		return {"R1":expand("FASTQ/TRIMMED/{sample}_{runID}_R1.fastq.gz",runID=runs, **wildcards),
		"R2":expand("FASTQ/TRIMMED/{sample}_{runID}_R2.fastq.gz",runID=runs, **wildcards)}



def star_input(wildcards):
	if not config["CUTADAPT"]["trimming"]:
		fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
		return {"R1": [f for f in fastqs.fastq1], "R2": [f for f in fastqs.fastq2]}
	else:
		fastqs = MASTER_LIST.loc[(wildcards.sample), ["fastq1", "fastq2"]].dropna()
		runs = fastqs.index.tolist()
		return {"R1":expand("FASTQ/TRIMMED/{sample}_{runID}_R1.fastq.gz",runID=runs, **wildcards),
		"R2":expand("FASTQ/TRIMMED/{sample}_{runID}_R2.fastq.gz",runID=runs, **wildcards)}