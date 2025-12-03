#!/bin/bash

vcftool_img=/mnt/pgx_src_data_pool_4/reference/images/consensus-call-anotate-v0.2.5.img
data_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/results/Mutect2
results_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/results/Mutect2
hc_region=/mnt/home_ssd/home/shangjun/project/fuscc_lc_1000/seqc2/wgs/variants_consistency/seqc2_files/High-Confidence_Regions_v1.2.bed

tnseq_dir=${data_dir}
singularity exec ${vcftool_img} vcftools --gzvcf ${tnseq_dir}/WES_IL_T_1.mutect2.filtered.vcf.gz --remove-indels --recode --out ${results_dir}/WES_IL_T_1.mutect2.filtered.snv
singularity exec ${vcftool_img} vcftools --gzvcf ${results_dir}/WES_IL_T_1.mutect2.filtered.snv.recode.vcf --bed ${hc_region}  --recode --out ${results_dir}/hc/WES_IL_T_1.mutect2.PASS.snv.hc
singularity exec ${vcftool_img} vcftools --gzvcf ${tnseq_dir}/WES_IL_T_1.mutect2.filtered.vcf.gz  --keep-only-indels --recode --out ${results_dir}/WES_IL_T_1.mutect2.filtered.indel
singularity exec ${vcftool_img} vcftools --gzvcf ${results_dir}/WES_IL_T_1.mutect2.filtered.indel.recode.vcf --bed ${hc_region}  --recode --out ${results_dir}/hc/WES_IL_T_1.mutect2.PASS.indel.hc