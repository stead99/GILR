#!/bin/bash

# images
# We only use the GATK image now, vcftool_img is deprecated in this script
# vcftool_img=/mnt/pgx_src_data_pool_4/reference/images/consensus-call-anotate-v0.2.5.img
IMG_GATK=/mnt/pgx_src_data_pool_4/reference/images/gatk-v4.6.2.0.img

# dir
data_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/results/Mutect2
results_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/results/Mutect2
hc_region=/mnt/home_ssd/home/shangjun/project/fuscc_lc_1000/seqc2/wgs/variants_consistency/seqc2_files/High-Confidence_Regions_v1.2.bed
tnseq_dir=${data_dir}
MUTECT2_FILTERED_VCF=${tnseq_dir}/WES_FD_T_1.mutect2.filtered.vcf.gz

# Ensure the output directory for 'hc' exists
mkdir -p "${results_dir}/hc"

# ----------------------------------------------------
# 1. Filter for SNPs and apply High-Confidence regions
# ----------------------------------------------------
echo "Filtering SNPs and applying high-confidence BED regions using GATK SelectVariants..."
# Output path: ${results_dir}/hc/WES_FD_T_1.mutect2.PASS.snv.hc.vcf.gz
singularity exec $IMG_GATK gatk SelectVariants \
     -V "${MUTECT2_FILTERED_VCF}" \
     --select-type-to-include SNP \
     -L "${hc_region}" \
     -O "${results_dir}/hc/WES_FD_T_1.mutect2.PASS.snv.hc.vcf.gz"


# ----------------------------------------------------
# 2. Filter for Indels (temporary file)
# ----------------------------------------------------
echo "Filtering for all indels using GATK SelectVariants..."
TEMP_INDEL_VCF="${results_dir}/WES_FD_T_1.mutect2.filtered.indel.vcf.gz"

singularity exec $IMG_GATK gatk SelectVariants \
     -V "${MUTECT2_FILTERED_VCF}" \
     --select-type-to-include INDEL \
     -O "${TEMP_INDEL_VCF}"


# ----------------------------------------------------
# 3. Apply High-Confidence regions to Indels
# ----------------------------------------------------
echo "Applying high-confidence BED regions to indels using GATK SelectVariants..."
# Output path: ${results_dir}/hc/WES_FD_T_1.mutect2.PASS.indel.hc.vcf.gz
singularity exec $IMG_GATK gatk SelectVariants \
     -V "${TEMP_INDEL_VCF}" \
     -L "${hc_region}" \
     -O "${results_dir}/hc/WES_FD_T_1.mutect2.PASS.indel.hc.vcf.gz"


# ----------------------------------------------------
# 4. Cleanup temporary files
# ----------------------------------------------------
echo "Cleaning up temporary indel VCF files..."
rm "${TEMP_INDEL_VCF}"
rm "${TEMP_INDEL_VCF}.tbi" # Remove the automatically generated index file
