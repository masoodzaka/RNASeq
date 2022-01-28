__author__ = "Masood Zaka (https://github.com/masoodzaka/RNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

configfile: "config.yaml"

# include functions.py modules
include: "SCRIPTS/functions.py"

STAR_INDEX = ["INDEX/STAR_INDEX"]
SALMON_INDEX = ["INDEX/SALMON_INDEX"]
STAR_BAMS= expand("STAR/{sample}Aligned.out.bam", sample=SAMPLES)
FEATURECOUNTS= expand("featureCounts/{sample}_featureCounts.txt", sample=SAMPLES)
HTSEQCOUNTS= expand("HTSeqCounts/{sample}_htseq.counts", sample=SAMPLES)
SALMONQUANTS= expand("SALMON/{sample}/quant.sf", sample=SAMPLES)
ARRIBA= expand("ARRIBA/{sample}.tsv", sample=SAMPLES)
SORTEDBAM=expand("SortedBAM/{sample}.sorted.bam", sample=SAMPLES)
SORTEDBAMINDEX=expand("SortedBAM/{sample}.sorted.bam.bai",sample=SAMPLES)
BIGWIG=expand("BIGWIG/{sample}.bw", sample=SAMPLES)



# extend all the rules using python function 
ALL = []

if config["STAR_INDEX"]:
    ALL.extend(STAR_INDEX)
ALL.extend(STAR_BAMS)

if config["SALMON_INDEX"]:
    ALL.extend(SALMON_INDEX)

if config["featureCounts"]:
    ALL.extend(FEATURECOUNTS)

if config["HTSeqCounts"]:
    ALL.extend(HTSEQCOUNTS)

#ALL.extend(SALMONQUANTS)

ALL.extend(ARRIBA)

ALL.extend(SORTEDBAM)

ALL.extend(SORTEDBAMINDEX)

ALL.extend(BIGWIG)

##

rule ALL:
    input:
        ALL,

rule Create_Star_Index:
    input:
        REF=config["REF"],
    
    output:
        INDEX=directory("INDEX/STAR_INDEX"),
    
    params:
        GTF=config["GTF"],

    log: "LOGS/STAR/star_index.log"
   
    benchmark:"LOGS/STAR/star_index.tsv"
    
    threads: 4
    
    conda: "ENVS/star.yaml"
    
    message:" Creating STAR index for {input} and saving as {output}"
    
    shell:""" STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.INDEX} \
            --genomeFastaFiles {input.REF} \
            --sjdbOverhang 100 \
            --sjdbGTFfile {params.GTF} 1&2>{log}
        
"""


rule Star_Aligner:
    input:
        unpack(star_input)

    output:
        BAM=("STAR/{sample}Aligned.out.bam"),
        TAB=("STAR/{sample}ReadsPerGene.out.tab")

    params:
        R1=lambda wildcards, input: str(",".join(input.R1)),
        R2=lambda wildcards, input: str(",".join(input.R2)),
        INDEX="INDEX/STAR_INDEX",
        PREFIX="STAR/{sample}"

    log: "LOGS/STAR/{sample}.log"

    benchmark:"LOGS/STAR/{sample}.tsv"

    threads: config["STAR_THREADS"]

    conda: "ENVS/star.yaml"

    message:" Running STAR aligner for {input} using {threads} threads and saving as {output}"

    shell:""" STAR --runMode alignReads \
    --runThreadN {threads} \
    --genomeDir {params.INDEX} \
    --genomeLoad NoSharedMemory \
    --readFilesIn {params.R1} {params.R2} \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --runRNGseed 777 \
    --outFilterType Normal \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --outFilterMultimapScoreRange 1 \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --outReadsUnmapped None \
    --alignIntronMin 20 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --sjdbOverhang 100 \
    --chimSegmentMin 20 \
    --chimJunctionOverhangMin 20 \
    --chimSegmentReadGapMax 3 \
    --chimOutType WithinBAM SoftClip \
    --quantMode GeneCounts \
    --outMultimapperOrder Random \
    --outSAMstrandField intronMotif \
    --outSAMattributes All \
    --outSAMunmapped Within KeepPairs \
    --outSAMtype BAM Unsorted \
    --limitBAMsortRAM 30000000000 \
    --outSAMmode Full \
    --outSAMheaderHD @HD VN:1.4 \
    --outFileNamePrefix {params.PREFIX} 2> {log}
"""

rule featureCounts:
    input:
        BAM=("STAR/{sample}Aligned.out.bam"),

    output:
        COUNTS=("featureCounts/{sample}_featureCounts.txt"),

    params:
       GTF=config["GTF"],

    log: "LOGS/featureCounts/{sample}.log"

    benchmark:"LOGS/featureCounts/{sample}.tsv"

    threads: 2

    conda: "ENVS/subread.yaml"

    message:" Running featureCounts for {input} using {threads} threads and saving as {output}"

    shell:"""
        featureCounts \
        -T {threads} \
        -a {params.GTF} \
        -p \
        -t exon \
        -g gene_id \
        -o {output.COUNTS} \
        {input.BAM} 2> {log}
"""

rule HTSeqCounts:
    input:
        BAM=("STAR/{sample}Aligned.out.bam"),

    output:
        COUNTS=("HTSeqCounts/{sample}_htseq.counts"),

    params:
       GTF=config["GTF"],

    log: "LOGS/HTSeqCounts/{sample}.log"

    benchmark:"LOGS/HTSeqCounts/{sample}.tsv"

    threads: 2

    conda: "ENVS/htseq.yaml"

    message:" Running HTSeqCounts for {input} using {threads} threads and saving as {output}"

    shell:"""
        htseq-count \
        -m intersection-nonempty \
        -i gene_id \
        -r name \
        -s no \
        -f bam \
        {input.BAM} \
        {params.GTF} > {output.COUNTS} 2> {log}
"""

rule Create_Salmon_Index:
    input:
        GTF=config["SALMON_GTF"],

    output:
        INDEX=directory("INDEX/SALMON_INDEX"),

    params:
        DECOYS=config["DECOYS"],
        extra="",

    log: "LOGS/SALMON/salmon_index.log"

    benchmark:"LOGS/SALMON/salmon_index.tsv"

    threads: 2

    conda: "ENVS/salmon.yaml"

    message:" Creating SALMON index for {input} and saving as {output}"

    shell:""" salmon \
            index \
            --threads {threads} \
            -t {input.GTF} \
            -d {params.DECOYS} \
            --gencode \
            -i {output.INDEX} 2>{log}

"""

rule Salmon_Quants:
    input:
        unpack(salmon_input)

    output:
        QUANT=("SALMON/{sample}/quant.sf"),
        LIB=("SALMON/{sample}/lib_format_counts.json")

    params:
        R1=lambda wildcards, input: str(" ".join(input.R1)),
        R2=lambda wildcards, input: str(" ".join(input.R2)),
        INDEX="INDEX/SALMON_INDEX",

    log: "LOGS/SALMON/{sample}.log"

    benchmark:"LOGS/SALMON/{sample}.tsv"

    threads: 2

    conda: "ENVS/salmon.yaml"

    message:" Running Salmon Quant for {input} using {threads} threads and saving as {output}"

    shell:""" salmon quant \
    -i {params.INDEX} \
    -l A \
    -1 {params.R1} \
    -2 {params.R2} \
    -p {threads} \
    --validateMappings \
    -o {output.QAUNT} 2> {log}
"""

rule Arriba:
    input:
        BAM=("STAR/{sample}Aligned.out.bam")

    output:
        FUSIONS=("ARRIBA/{sample}.tsv"),
        DISCARDED=("ARRIBA/{sample}.discarded.tsv")

    params:
        GTF=config["GTF"],
        REF=config["REF"],
        BLACKLIST=config["BLACKLIST"],
        KNOWNFUSION=config["KNOWNFUSION"],

    log: "LOGS/ARRIBA/{sample}.log"

    benchmark:"LOGS/ARRIBA/{sample}.tsv"

    threads: 2

    conda: "ENVS/arriba.yaml"

    message:" Running Arriba Fusion for {input} using {threads} threads and saving as {output}"

    shell:""" arriba \
    -x {input.BAM} \
    -a {params.REF} \
    -g {params.GTF} \
    -b {params.BLACKLIST} \
    -k {params.KNOWNFUSION} \
    -t {params.KNOWNFUSION} \
    -o {output.FUSIONS} \
    -O {output.DISCARDED} 2> {log}
"""
rule SortedBAM:
    input:
        BAM=("STAR/{sample}Aligned.out.bam"),

    output:
        BAM=("SortedBAM/{sample}.sorted.bam"),

    params:
        MEM="2G",

    log: "LOGS/SortedBAM/{sample}.log"

    benchmark: "LOGS/SortedBAM/{sample}.tsv"

    threads: 2

    conda:"ENVS/arriba.yaml"

    message: "Running Samtools Sort for {input} using {threads} threads and saving as {output}"

    shell:""" samtools \
    sort \
    -m {params.MEM} \
    -@ 5 \
    -T {output.BAM} \
    -o {output.BAM} \
    {input.BAM} 2> {log}
    """

rule SortedBAMIndex:
    input:
        BAM=("SortedBAM/{sample}.sorted.bam"),

    output:
        BAI=("SortedBAM/{sample}.sorted.bam.bai")

    params:
        MEM="2G",

    log: "LOGS/SortedBAM/{sample}.index.log"

    benchmark: "LOGS/SortedBAM/{sample}.index.tsv"

    threads: 2

    conda:"ENVS/arriba.yaml"

    message: "Running Samtools Index for {input} using {threads} threads and saving as {output}"

    shell:""" samtools \
    index \
    {input.BAM} 2> {log}
    """

rule CreateBigWig:
    input:
        BAM=("SortedBAM/{sample}.sorted.bam"),
        BAI=("SortedBAM/{sample}.sorted.bam.bai"),

    output:
        BW=("BIGWIG/{sample}.bw")

    params:
        BINSIZE=20,

    log: "LOGS/BIGWIG/{sample}.log"

    benchmark: "LOGS/BIGWIG/{sample}.tsv"

    threads: 2

    conda:"ENVS/deeptools.yaml"

    message: "Running Deeptools bamCoverage for {input} using {threads} threads and saving as {output}"

    shell:""" bamCoverage \
    -b {input.BAM} \
    -p {threads} \
    -bs {params.BINSIZE} \
    --smoothLength 60 \
    --skipNonCoveredRegions \
    --normalizeUsing RPKM \
    -o {output.BW} 2> {log}
    """