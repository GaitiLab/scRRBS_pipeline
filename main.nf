#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define workflow input parameters
params.reads = "{sample}/fastq/{pool}_L00*_{1,2}_001.fastq.gz"
params.scripts = './bin' // Default path, will be overridden by config

workflow {

    // Create a channel for input files and group by sample and pool
    Channel
        .fromFilePairs(params.reads, size: 2, flat: true)
        .map { sample, files -> tuple(sample.split('/')[0], sample.split('/')[2].split('_')[0], files) }
        .set { fastq_files }

    // Process to merge lanes
    process mergeLanes {
        input:
        tuple val(sample), val(pool), path(reads)

        output:
        path("{sample}/merged/{pool}_ME_L001_R1_001.fastq.gz"), emit: mergedRead1
        path("{sample}/merged/{pool}_ME_L001_R2_001.fastq.gz"), emit: mergedRead2

        script:
        """
        ${params.scripts}/merge_lanes.sh {sample}/merged/{pool}_ME_L001_R1_001.fastq.gz ${reads[0]}
        ${params.scripts}/merge_lanes.sh {sample}/merged/{pool}_ME_L001_R2_001.fastq.gz ${reads[1]}
        """
    }

    // Execute the process
    mergeLanes(fastq_files)
}
