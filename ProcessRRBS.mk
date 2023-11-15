import os
import glob

configfile: "config.yaml"

sample_dirs = [d for d in next(os.walk('.'))[1] if d != "logs" and not d.startswith('.') and d != "slurm_logs"]

# Function to get cells for each sample
def get_pools(sample):
    fastq_files = glob.glob(os.path.join(sample, "fastq", "*.fastq.gz"))
    pools = set(os.path.basename(f).split("_L00")[0] for f in fastq_files)
    return list(pools)

# Collect all expected output files for the `all` rule
merged_files = []
for sample in sample_dirs:
    for pool in get_pools(sample):
        for read in [1, 2]:
            merged_files.append(f"{sample}/merged/{pool}_ME_L001_R{read}_001.fastq")

demultiplexed_files = []
for sample in sample_dirs:
	for pool in get_pools(sample):
		for barcode in config["barcodes"]:
			for read in [1, 2]:
				demultiplexed_files.append(f"{sample}/demultiplexed/{pool}.{barcode}.R{read}.fastq")

aligned_files = []
for sample in sample_dirs:
	for pool in get_pools(sample):
		for barcode in config["barcodes"]:
			aligned_files.append(f"{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_pe.bam")

extracted_files = []
for sample in sample_dirs:
	for pool in get_pools(sample):
		for barcode in config["barcodes"]:
			extracted_files.append(f"{sample}/extracted/{pool}.{barcode}.R1_bismark_bt2_pe.bismark.cov.gz")
			extracted_files.append(f"{sample}/extracted/CpG_context_{pool}.{barcode}.R1_bismark_bt2_pe.txt")

rule all:
	input:
		merged_files,
		demultiplexed_files,
		aligned_files

# rule sites_only:
# 	input:
# 		"site_matrix.csv",
#         "Pairwise_stats_rrbs.csv"

# rule by_region:
# 	input:
# 		"region_matrix_pct.csv",
# 		"region_matrix_count.csv",
# 		"Pairwise_stats_rrbs.csv"

rule merge_lanes:
    input:
        r1=lambda wildcards: sorted(glob.glob(f"{wildcards.sample}/fastq/{wildcards.pool}_L00*_R1_001.fastq.gz")),
        r2=lambda wildcards: sorted(glob.glob(f"{wildcards.sample}/fastq/{wildcards.pool}_L00*_R2_001.fastq.gz"))
    output:
        r1="{sample}/merged/{pool}_ME_L001_R1_001.fastq",
        r2="{sample}/merged/{pool}_ME_L001_R2_001.fastq"
    log:
        "logs/merge_lanes_{sample}_{pool}.log"
    shell:
        """
        {config[scripts]}/merge_lanes.sh {output.r1} {input.r1}
        {config[scripts]}/merge_lanes.sh {output.r2} {input.r2}
        """

rule demultiplex:
	input:
		r1="{sample}/merged/{pool}_ME_L001_R1_001.fastq",
		r2="{sample}/merged/{pool}_ME_L001_R2_001.fastq"
	output:
		expand("{{sample}}/demultiplexed/{{pool}}.{barcode}.R1.fastq", barcode=config["barcodes"]),
		expand("{{sample}}/demultiplexed/{{pool}}.{barcode}.R2.fastq", barcode=config["barcodes"])
	shell:
		"""
		output_prefix=`python -c "print('{input.r1}'.replace('_ME_L001_R1_001.fastq', '').replace('merged','demultiplexed'))"`
		perl {config[scripts]}/splitFastqPair.pl {input.r1} {input.r2} $output_prefix
		"""

rule bismark_align:
	input:
		r1="{sample}/demultiplexed/{pool}.{barcode}.R1.fastq",
		r2="{sample}/demultiplexed/{pool}.{barcode}.R2.fastq"
	output:
		"{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_pe.bam"
	threads: 4
	shell:
		"bismark --multicore {threads} -X 1000 --path_to_bowtie $(dirname $(which bowtie2)) --un --ambiguous {config[genome]} -1 {input.r1} -2 {input.r2} -o {wildcards.sample}/aligned/ --temp_dir {wildcards.sample}/aligned/"

rule bismark_extract:
	input:
		"{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_pe.bam"
	output:
		"{sample}/extracted/{pool}.{barcode}.R1_bismark_bt2_pe.bismark.cov.gz",
		"{sample}/extracted/CpG_context_{pool}.{barcode}.R1_bismark_bt2_pe.txt"
	threads: 4
	shell:
		"bismark_methylation_extractor --multicore {threads} --bedgraph --comprehensive {input} -o {wildcards.sample}/extracted/"

# rule pdr:
# 	input:
# 		expand("CpG_context_{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.txt", prefix=config["prefixes"], barcode=config["barcodes"])
# 	output:
# 		"PDR_rrbs.txt"
# 	shell:
# 		"python {config[scripts]}/compute_pdr.py {input} > {output}"

# rule qc_stats:
# 	input:
# 		"PDR_rrbs.txt"
# 	params:
# 		cells=expand("{prefix}.{barcode}", prefix=config["prefixes"], barcode=config["barcodes"])
# 	output:
# 		"QC_stats_rrbs.csv"
# 	shell:
# 		"python {config[scripts]}/cov_batch.py {params.cells} {input} rrbs"

# rule pairwise_stats:
# 	input:
# 		"QC_stats_rrbs.csv"
# 	output:
# 		"Pairwise_stats_rrbs.csv"
# 	shell:
# 		"python {config[scripts]}/compute_pair.py {input} > {output}"

# rule site_matrix:
# 	input:
# 		expand("{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.bismark.cov.gz", prefix=config["prefixes"], barcode=config["barcodes"])
# 	output:
# 		"site_matrix.csv"
# 	shell:
# 		"python {config[scripts]}/build_features.py {input} > {output}"

# rule site_bed:
# 	input:
# 		"site_matrix.csv"
# 	output:
# 		"site_matrix.bed"
# 	shell:
# 		"python {config[scripts]}/csv_to_site_bed.py {input} > {output}"

# rule site_intersect:
# 	input:
# 		"site_matrix.bed"
# 	output:
# 		"site_intersect.bed"
# 	shell:
# 		"intersectBed -a {config[regions]} -b {input} -wa -wb > {output}"

# rule region_matrix:
# 	input:
# 		"site_matrix.csv",
# 		"site_intersect.bed"
# 	output:
# 		"region_matrix_pct.csv",
# 		"region_matrix_count.csv"
# 	shell:
# 		"python {config[scripts]}/met_matrix.py {config[path]}"
