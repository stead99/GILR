#!/bin/bash

# ------------------------------------------------
# Function to process a single sample
# Argument 1: Sample Name (e.g., U1a_ACAGTG_L001)
# 
# IMPORTANT: This function relies on global variables 
# defined and sourced from the run_pipeline.sh script:
# data_dir, results_base, log_dir, fasta, R1_suffix, R2_suffix, 
# IMG_BWA, IMG_SAMTOOLS, IMG_GATK, nt, pl.
# ------------------------------------------------
process_sample() {
    nam=$1
    # Set up sample-specific log file and redirect stderr
    logfile="$nam.$(date +"%Y%m%d%H%M").log"
    log_message() {
        echo "$1" | tee -a "$logfile"
    }
    
    # 确保日志文件被创建（清空旧内容，如果有的话）
    > "$logfile"
    set -x 2>&1 | tee -a "$logfile" # 将 set -x 的输出（命令本身）也捕获到日志

    log_message "--- Starting processing for sample: $nam ---" 

    # Define input file paths
    fastq_1=${data_dir}/${nam}${R1_suffix}
    fastq_2=${data_dir}/${nam}${R2_suffix}
    group=$nam # Read group ID
    sample=$nam # Sample name

    # Define module-specific output directories for this sample
    mapping_dir=${results_base}/mapping
    dedup_dir=${results_base}/dedup
    mkdir -p "$mapping_dir" "$dedup_dir"

    # ******************************************
    # 1. Mapping reads with BWA-MEM & Sorting with Samtools
    # ******************************************
    start_time_mod=$(date +%s)
    singularity exec $IMG_BWA bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$pl" \
        -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 | \
        singularity exec $IMG_SAMTOOLS samtools sort -o ${mapping_dir}/${nam}.sorted.bam
    
    # build index
    singularity exec $IMG_SAMTOOLS samtools index ${mapping_dir}/${nam}.sorted.bam

    end_time_mod=$(date +%s)
    echo "Module BWA+Sort Elapsed time: "$(($end_time_mod - $start_time_mod))" sec" >>${log_dir}/$logfile

    # ******************************************
    # 2. Remove Duplicate Reads with GATK MarkDuplicates
    # ******************************************
    start_time_mod=$(date +%s)
    singularity exec $IMG_GATK gatk MarkDuplicates \
      -I ${mapping_dir}/${nam}.sorted.bam \
      -O ${dedup_dir}/${nam}.dedup.bam \
      -M ${dedup_dir}/${nam}.duplication.metrics \
      --REMOVE_DUPLICATES true \
      --CREATE_INDEX true

    end_time_mod=$(date +%s)
    echo "Module dedup Elapsed time: "$(($end_time_mod - $start_time_mod))" sec" >>${log_dir}/$logfile
    echo "--- Finished processing for sample: $nam ---" >&2
}

# Export the function so it is available to the main script when sourced
export -f process_sample