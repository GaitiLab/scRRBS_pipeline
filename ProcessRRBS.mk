import os
import glob

configfile: "config.yaml"

# sample_dirs = [d for d in next(os.walk('.'))[1] if d != "logs" and not d.startswith('.') and d != "slurm_logs"]

# # Function to get cells for each sample
# def get_pools(sample):
#     fastq_files = glob.glob(os.path.join(sample, "fastq", "*.fastq.gz"))
#     pools = set(os.path.basename(f).split("_L00")[0] for f in fastq_files)
#     return list(pools)

# # Collect all expected output files for the `all` rule
# merged_files = []
# for sample in sample_dirs:
#     for pool in get_pools(sample):
#         for read in [1, 2]:
#             merged_files.append(f"{sample}/merged/{pool}_ME_L001_R{read}_001.fastq")

# demultiplexed_files = []
# for sample in sample_dirs:
#     for pool in get_pools(sample):
#         for barcode in config["barcodes"]:
#             for read in [1, 2]:
#                 demultiplexed_files.append(f"{sample}/demultiplexed/{pool}.{barcode}.R{read}.fastq")

# aligned_files = []
# for sample in sample_dirs:
#     for pool in get_pools(sample):
#         for barcode in config["barcodes"]:
#             aligned_files.append(f"{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_pe.bam")

# status_files = []
# for sample in sample_dirs:
#     for pool in get_pools(sample):
#         for barcode in config["barcodes"]:
#             status_files.append(f"{sample}/aligned/{pool}.{barcode}.efficiency_checked")

# rule all:
#     input:
#         merged_files,
#         demultiplexed_files,
#         aligned_files,
#         status_files,
#         expand("{sample}/site_matrix.csv", sample=sample_dirs),
#         expand("{sample}/PDR_rrbs.txt", sample=sample_dirs)        

# rule sites_only:
#     input:
#         "site_matrix.csv",
#         "Pairwise_stats_rrbs.csv"

# rule by_region:
#     input:
#         "region_matrix_pct.csv",
#         "region_matrix_count.csv",
#         "Pairwise_stats_rrbs.csv"

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
    log:
        "logs/demultiplex_{sample}_{pool}.log"
    shell:
        """
        output_prefix=`python -c "print('{input.r1}'.replace('_ME_L001_R1_001.fastq', '').replace('merged','demultiplexed'))"`
        perl {config[scripts]}/splitFastqPair.pl {input.r1} {input.r2} $output_prefix
        """

# rule bismark_align:
#     input:
#         r1="{sample}/demultiplexed/{pool}.{barcode}.R1.fastq",
#         r2="{sample}/demultiplexed/{pool}.{barcode}.R2.fastq"
#     output:
#         "{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_pe.bam",
#         "{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_PE_report.txt"
#     log:
#         "logs/bismark_align_{sample}_{pool}_{barcode}.log"
#     threads: 4
#     shell:
#         "bismark --multicore {threads} -X 1000 --path_to_bowtie $(dirname $(which bowtie2)) --un --ambiguous {config[genome]} -1 {input.r1} -2 {input.r2} -o {wildcards.sample}/aligned/ --temp_dir {wildcards.sample}/aligned/"

# rule check_alignment_efficiency:
#     input:
#         report="{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_PE_report.txt"
#     output:
#         marker=temp("{sample}/aligned/{pool}.{barcode}.efficiency_checked")
#     params:
#         cutoff=config['efficiency']
#     log:
#         "logs/check_alignment_efficiency_{sample}_{pool}_{barcode}.log"
#     run:
#         with open(input.report, 'r') as report_file:
#             for line in report_file:
#                 if 'Mapping efficiency:' in line:
#                     efficiency = float(line.split(':')[1].strip().strip('%'))
#                     break
#         status = 'ok' if efficiency >= params.cutoff else 'no'
#         with open(f"{wildcards.sample}/aligned/{wildcards.pool}.{wildcards.barcode}.{status}", 'w') as status_file:
#             status_file.write('')
#         with open(output.marker, 'w') as marker_file:
#             marker_file.write('')

# checkpoint bismark_extract:
#     input:
#         ok="{sample}/aligned/{pool}.{barcode}.ok",
#         bam="{sample}/aligned/{pool}.{barcode}.R1_bismark_bt2_pe.bam"
#     output:
#         "{sample}/extracted/{pool}.{barcode}.R1_bismark_bt2_pe.bismark.cov.gz",
#         "{sample}/extracted/CpG_context_{pool}.{barcode}.R1_bismark_bt2_pe.txt"
#     log:
#         "logs/bismark_extract_{sample}_{pool}_{barcode}.log"
#     shell:
#         "bismark_methylation_extractor --bedgraph --comprehensive {input.bam} -o {wildcards.sample}/extracted/"

# from snakemake.io import expand

# def get_cov_files(wildcards):
#     ok_files = glob.glob("*/*/*.ok")
#     cov_files = []
#     for ok_file in ok_files:
#         parts = ok_file.split('/')
#         sample = parts[0]
#         pool, barcode = parts[2].split('.')[0:2]
#         cov_file = f"{sample}/extracted/{pool}.{barcode}.R1_bismark_bt2_pe.bismark.cov.gz"
#         cov_files.append(cov_file)
#     return cov_files

# # Use the function in the rule
# rule site_matrix:
#     input:
#         cov=lambda wildcards: sorted(glob.glob(f"{wildcards.sample}/extracted/{wildcards.pool}.{wildcards.barcode}.R1_bismark_bt2_pe.bismark.cov.gz")),
#         ok=lambda wildcards: sorted(glob.glob(f"{wildcards.sample}/aligned/{wildcards.pool}.{wildcards.barcode}.ok"))
#     output:
#         "{sample}/site_matrix.csv"
#     log:
#         "logs/site_matrix_{sample}.log"
#     shell:
#         "python {config[scripts]}/build_features.py {wildcards.sample} {input.cov}"

# def get_extracted_txt_files(wildcards):
#     ok_files = glob.glob("*/*/*.ok")
#     txt_files = []
#     for ok_file in ok_files:
#         parts = ok_file.split('/')
#         sample = parts[0]
#         pool, barcode = parts[2].split('.')[0:2]
#         txt_file = f"{sample}/extracted/CpG_context_{pool}.{barcode}.R1_bismark_bt2_pe.txt"
#         txt_files.append(txt_file)
#     return txt_files

# rule pdr:
#     input:
#         extracted_txt=get_extracted_txt_files
#     output:
#         "{sample}/PDR_rrbs.txt"
#     log:
#         "logs/pdr_{sample}.log"
#     shell:
#         "python {config[scripts]}/compute_pdr.py {wildcards.sample} {input}"

# rule qc_stats:
#     input:
#         "PDR_rrbs.txt"
#     params:
#         cells=expand("{prefix}.{barcode}", prefix=config["prefixes"], barcode=config["barcodes"])
#     output:
#         "QC_stats_rrbs.csv"
#     shell:
#         "python {config[scripts]}/cov_batch.py {params.cells} {input} rrbs"

# rule pairwise_stats:
#     input:
#         "QC_stats_rrbs.csv"
#     output:
#         "Pairwise_stats_rrbs.csv"
#     shell:
#         "python {config[scripts]}/compute_pair.py {input} > {output}"



# rule site_bed:
#     input:
#         "site_matrix.csv"
#     output:
#         "site_matrix.bed"
#     shell:
#         "python {config[scripts]}/csv_to_site_bed.py {input} > {output}"

# rule site_intersect:
#     input:
#         "site_matrix.bed"
#     output:
#         "site_intersect.bed"
#     shell:
#         "intersectBed -a {config[regions]} -b {input} -wa -wb > {output}"

# rule region_matrix:
#     input:
#         "site_matrix.csv",
#         "site_intersect.bed"
#     output:
#         "region_matrix_pct.csv",
#         "region_matrix_count.csv"
#     shell:
#         "python {config[scripts]}/met_matrix.py {config[path]}"
