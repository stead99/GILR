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
    
    # Check if 'tee' is available to decide logging method
    if command -v tee >/dev/null 2>&1; then
        # If tee exists, use a simple trick to pipe output, still might not work with extremely old shell
        echo "Using tee for logging to $logfile"
        # This part is tricky to do inside a function without process substitution
        # We fall back to redirecting ONLY to the file below if process substitution fails
    else
        echo "tee command not found or process substitution unsupported. Logging only to file: $logfile"
    fi

    # Redirect all subsequent output of this function ONLY to the log file (no screen output)
    # This is a highly compatible method but lacks simultaneous console output.
    exec > "$logfile" 2>&1

    set -x # Log commands as they are executed (these logs now only go to the file)


    # Define paths specific to this pair
    dedup_dir=${seqc_bam_dir} # Use the global data_dir as input source
    BQSR_dir=${results_base}/BQSR
    Mutect2_dir=${results_base}/Mutect2
    mkdir -p "$BQSR_dir" "$Mutect2_dir"

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

    # Tumor BQSR
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
    # Note: Placeholder 'pon.vcf.gz' and 'common_biallelic.vcf.gz' must exist in $gatk_res
    start_time_mod=$(date +%s)
    singularity exec $IMG_GATK gatk Mutect2 \
         -R "${fasta}" \
         -I "${TUMOR_RECAL_BAM}" \
         -I "${NORMAL_RECAL_BAM}" \
         -normal "$normal_name" \
         --germline-resource "${germline_resource}" \
         #--panel-of-normals "${pon_dir}/pon.vcf.gz" \
         -O "${MUTECT2_VCF}"

    ### GetPileupSummaries & CalculateContamination & FilterMutectCalls
    singularity exec $IMG_GATK gatk GetPileupSummaries \
       -I "${TUMOR_RECAL_BAM}" \
       --intervals "${bed_file}" \
       -V "${common_biallelic}" \
       -L "${common_biallelic}" \
       -O ${Mutect2_dir}/${tumor_name}.pileups.table

    singularity exec $IMG_GATK gatk GetPileupSummaries \
       -I "${NORMAL_RECAL_BAM}" \
       --intervals "${bed_file}" \
       -V "${common_biallelic}" \
       -L "${common_biallelic}" \
       -O ${Mutect2_dir}/${normal_name}.pileups.table

    singularity exec $IMG_GATK gatk CalculateContamination \
       -I ${Mutect2_dir}/${tumor_name}.pileups.table \
       -matched ${Mutect2_dir}/${normal_name}.pileups.table \
       -O ${Mutect2_dir}/${tumor_name}.contamination.table \
       --segments ${Mutect2_dir}/${tumor_name}.segments.tsv

    singularity exec $IMG_GATK gatk FilterMutectCalls \
       -R "${fasta}" \
       -V "${MUTECT2_VCF}" \
       --contamination-table ${Mutect2_dir}/${tumor_name}.contamination.table \
       --tumor-segmentation ${Mutect2_dir}/${tumor_name}.segments.tsv \
       -O ${Mutect2_dir}/${tumor_name}.mutect2.filtered.vcf.gz

    end_time_mod=$(date +%s)
    echo "Module Mutect2 Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"
    echo "--- Finished processing for pair: $pair_id ---"
    
    set +x
}

# Export the function so it is available to the main script when sourced
export -f process_tn_pair
