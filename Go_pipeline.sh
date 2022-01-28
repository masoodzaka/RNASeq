#!/bin/bash
set -e
set -u
set -o pipefail

CORES=$1

echo -e "\n"

echo "Masood Zaka 2021-2022"

echo -e "\nSubmitting Snakemake for RNA Analysis using $CORES cores\n"

echo ""

snakemake \
--use-conda \
--rerun-incomplete \
--cores $CORES \
--latency-wait 120 \
--printshellcmds

echo -e "\nAll jobs completed\n"