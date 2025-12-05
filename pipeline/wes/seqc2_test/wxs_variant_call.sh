#!/bin/bash

# ------------------------------------------------
# Function to process a single Tumor-Normal pair
# Argument 1: Tumor Sample Name prefix (e.g., WES_IL_T_1.bwa)
# Argument 2: Normal Sample Name prefix (e.g., WES_IL_N_1.bwa)
# 
# IMPORTANT: This function relies on global variables 
# defined and sourced from the run_tn_pipeline.sh script.
# ------------------------------------------------
process_tn_pair() {
    tumor_name=$1
    normal_name=$2
    pair_id="${tumor_name}_vs_${normal_name}" # Used for logging

    # Set up sample-specific log file (basic logging using simple redirection)
    logfile="${log_dir}/${pair_id}.$(date +"%Y%m%d%H%M").log"
    
    # --------------------------------------------------
    # 兼容性日志记录方法 (Compatibility Logging Method)
    # 将所有输出重定向到日志文件，不支持同时输出到屏幕
    # Redirects all output to the log file ONLY. 
    # Does not support simultaneous console output.
    # --------------------------------------------------
    echo "Logging all output for pair $pair_id to $logfile (console output disabled for compatibility)."
    exec > "$logfile" 2>&1

    set -x # Log commands as they are executed

    # Define paths specific to this pair
    dedup_dir=${results_base}/dedup
    BQSR_dir=${results_base}/BQSR
    Mutect2_dir=${results_base}/Mutect2
    strelka_dir=${results_base}/strelka
    mkdir -p "$BQSR_dir" "$Mutect2_dir" "$strelka_dir"

    TUMOR_BAM=${dedup_dir}/${tumor_name}.dedup.bam
    NORMAL_BAM=${dedup_dir}/${normal_name}.dedup.bam

    # Check if input BAMs exist
    if [[ ! -f "$TUMOR_BAM" || ! -f "$NORMAL_BAM" ]]; then
        echo "Error: Missing input BAM files for $pair_id. Skipping."
        set +x
        return 1
    fi

    # ******************************************
    # 1. BQSR (Base Quality Score Recalibration)
    # ******************************************
    start_time_mod=$(date +%s)
    
    --- Dynamic BQSR Commands (if you eventually uncomment them) ---
    These currently rely on fixed $db_mills and $dbsnp variables

    Tumor BQSR
    singularity exec $IMG_GATK gatk BaseRecalibrator \
      -I "$TUMOR_BAM" \
      -R "$fasta" \
      --known-sites "$db_mills" \
      --known-sites "$dbsnp" \
      -O ${BQSR_dir}/${tumor_name}.bqsr.grp

    singularity exec $IMG_GATK gatk ApplyBQSR \
      -R "$fasta" \
      -I "$TUMOR_BAM" \
      --bqsr-recal-file ${BQSR_dir}/${tumor_name}.bqsr.grp \
      -O ${BQSR_dir}/${tumor_name}.recalibrated.bam

    # Normal BQSR
    singularity exec $IMG_GATK gatk BaseRecalibrator \
      -I "$NORMAL_BAM" \
      -R "$fasta" \
      --known-sites "$db_mills" \
      --known-sites "$dbsnp" \
      -O ${BQSR_dir}/${normal_name}.bqsr.grp

    singularity exec $IMG_GATK gatk ApplyBQSR \
      -R "$fasta" \
      -I "$NORMAL_BAM" \
      --bqsr-recal-file ${BQSR_dir}/${normal_name}.bqsr.grp \
      -O ${BQSR_dir}/${normal_name}.recalibrated.bam

    end_time_mod=$(date +%s)
    echo "Module BQSR Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"

    # Define recalibrated BAM paths for subsequent steps
    TUMOR_RECAL_BAM=${BQSR_dir}/${tumor_name}.recalibrated.bam
    NORMAL_RECAL_BAM=${BQSR_dir}/${normal_name}.recalibrated.bam
    MUTECT2_VCF=${Mutect2_dir}/${tumor_name}.mutect2.unfiltered.vcf.gz

    # ******************************************
    # 2. Mutect2 Somatic Variant Calling
    # ******************************************
    start_time_mod=$(date +%s)
    
    # --- Dynamic Mutect2 Command Building ---
    MUTECT2_CMD="singularity exec $IMG_GATK gatk Mutect2 -R \"${fasta}\" -I \"${TUMOR_RECAL_BAM}\" -I \"${NORMAL_RECAL_BAM}\""
    
    if [[ -f "${bed_file}" ]]; then
        echo "Adding --intervals flag to Mutect2 using ${bed_file}"
        MUTECT2_CMD+=" --intervals \"${bed_file}\""
    else
        echo "Warning: bed_file not found. Mutect2 running whole-genome."
    fi

    MUTECT2_CMD+=" -normal \"$normal_name\" --germline-resource \"${germline_resource}\" -O \"${MUTECT2_VCF}\""
    eval $MUTECT2_CMD
    # ----------------------------------------

    ### GetPileupSummaries & CalculateContamination & FilterMutectCalls

    # --- Dynamic GetPileupSummaries Command Building ---
    # Use dynamic command building to omit the interval flags completely if the file is missing
    GPS_TUMOR_CMD="singularity exec $IMG_GATK gatk GetPileupSummaries -I \"${TUMOR_RECAL_BAM}\" -V \"${common_biallelic}\" -L \"${common_biallelic}\" -O ${Mutect2_dir}/${tumor_name}.pileups.table"
    GPS_NORMAL_CMD="singularity exec $IMG_GATK gatk GetPileupSummaries -I \"${NORMAL_RECAL_BAM}\" -V \"${common_biallelic}\" -L \"${common_biallelic}\" -O ${Mutect2_dir}/${normal_name}.pileups.table"
    
    if [[ -f "${bed_file}" ]]; then
        echo "Adding --intervals and -L flags to GetPileupSummaries using ${bed_file}"
        GPS_TUMOR_CMD+=" --intervals \"${bed_file}\""
        GPS_NORMAL_CMD+=" --intervals \"${bed_file}\""
    else
        echo "Warning: bed_file not found. GetPileupSummaries running whole-genome (may take longer)."
    fi
    
    echo "Executing Tumor GetPileupSummaries command..."
    eval $GPS_TUMOR_CMD

    echo "Executing Normal GetPileupSummaries command..."
    eval $GPS_NORMAL_CMD
    # --------------------------------------------------

    echo "Calculating Contamination..."
    singularity exec $IMG_GATK gatk CalculateContamination \
       -I ${Mutect2_dir}/${tumor_name}.pileups.table \
       -matched ${Mutect2_dir}/${normal_name}.pileups.table \
       -O ${Mutect2_dir}/${tumor_name}.contamination.table \
       --segments ${Mutect2_dir}/${tumor_name}.segments.tsv

    echo "Filtering Mutect Calls..."
    singularity exec $IMG_GATK gatk FilterMutectCalls \
       -R "${fasta}" \
       -V "${MUTECT2_VCF}" \
       --contamination-table ${Mutect2_dir}/${tumor_name}.contamination.table \
       --tumor-segmentation ${Mutect2_dir}/${tumor_name}.segments.tsv \
       -O ${Mutect2_dir}/${tumor_name}.mutect2.filtered.vcf.gz
    
    # ******************************************
    # 4. Strelka
    # ******************************************
    start_time_mod=$(date +%s)

    # --- Dynamic Strelka Configure Command Building ---
    STRELKA_CFG_CMD="singularity exec $IMG_strelka configureStrelkaSomaticWorkflow.py"
    STRELKA_CFG_CMD+=" --referenceFasta ${fasta} --tumorBam ${TUMOR_RECAL_BAM} --normalBam ${NORMAL_RECAL_BAM} --runDir ${strelka_dir}/${pair_id}"
    
    if [[ -f "${ex_decoys}" ]]; then
        echo "Adding --callRegions flag to Strelka using ${ex_decoys}"
        STRELKA_CFG_CMD+=" --callRegions ${ex_decoys}"
    else
        echo "Warning: ex_decoys not found. Strelka running whole-genome."
    fi
    
    eval $STRELKA_CFG_CMD
    # --------------------------------------------------

    singularity exec $IMG_strelka \
    ${strelka_dir}/${pair_id}/runWorkflow.py -m local -j $nt

    zcat ${strelka_dir}/${pair_id}/results/variants/somatic.snvs.vcf.gz|awk -F '\t' '{if(($1~"^#")||($1!~"^#" && $7=="PASS")){print $0}}' > ${strelka_dir}/${pair_id}/results/variants/${pair_id}.snvs.PASS.vcf
    zcat ${strelka_dir}/${pair_id}/results/variants/somatic.indels.vcf.gz|awk -F '\t' '{if(($1~"^#")||($1!~"^#" && $7=="PASS")){print $0}}' > ${strelka_dir}/${pair_id}/results/variants/${pair_id}.indels.PASS.vcf

    end_time_mod=$(date +%s)
    echo "Module strelka Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"
    echo "--- Finished processing for pair: $pair_id ---"
    set +x
}

# Export the function so it is available to the main script when sourced
export -f process_tn_pair
