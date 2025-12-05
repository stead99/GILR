#!/bin/bash

# ------------------------------------------------
# 全局配置变量 (Global Configuration Variables)
# ------------------------------------------------

# Base input data directory (where .dedup.bam files are)
data_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq

# Base output directory for results
results_base=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/results
pon_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/results/pon

# Log directory
log_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/log

# Reference paths
fasta=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/gatk_resource/Homo_sapiens_assembly38.fasta
gatk_res=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/gatk_resource
germline_resource=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/sentieon/germline_resource/af_info_gnomad_genomes_v3.1.2.vcf.gz
dbsnp=$gatk_res/Homo_sapiens_assembly38.dbsnp138.vcf
db_mills=$gatk_res/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
common_biallelic=$gatk_res/somatic-hg38_small_exac_common_3.hg38.vcf.gz
bed_file=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/bed_files/novogene_agilent_region.hg38.bed
ex_decoys=//mnt/pgx_src_data_pool_4/reference/human/GRCh38/bed_files/novogene_agilent_region.hg38.bed.gz
# Note: common_biallelic.vcf.gz and pon.vcf.gz must also exist in $gatk_res

# Execution environments (Singularity images)
IMG_GATK=/mnt/pgx_src_data_pool_4/reference/images/gatk-v4.6.2.0.img

# Number of threads (nt is less critical for GATK in single-sample mode but kept for consistency)
nt=16

# Ensure log directory exists
mkdir -p "$log_dir"

# ------------------------------------------------
# 引入核心处理函数 (Source the core processing function)
# ------------------------------------------------
source ./wxs_variant_call.sh

# ------------------------------------------------
# 从外部文件读取样本对列表 (Read sample pairs list from an external file)
# ------------------------------------------------

# 定义样本对列表文件的路径和名称
SAMPLE_PAIRS_FILE="./sample_pairs_il.txt"

# 检查样本列表文件是否存在
if [[ ! -f "$SAMPLE_PAIRS_FILE" ]]; then
    echo "错误：找不到样本列表文件 '$SAMPLE_PAIRS_FILE'。"
    echo "请创建一个名为 sample_pairs.txt 的文件，并在其中每行放置一对样本名称 (Tumor Normal，用空格或制表符分隔)。"
    exit 1
fi

echo "正在从 $SAMPLE_PAIRS_FILE 读取样本对列表..."

# ------------------------------------------------
# 循环处理所有样本对 (Loop through all sample pairs)
# ------------------------------------------------

# 使用 while read 循环读取文件内容
# IFS 设置为空格和制表符，以便正确解析 TUMOR_PREFIX 和 NORMAL_PREFIX
while IFS=$' \t' read -r TUMOR_PREFIX NORMAL_PREFIX || [[ -n "$TUMOR_PREFIX" ]]; do
    # 跳过空行或注释行（以 # 开头）
    if [[ -z "$TUMOR_PREFIX" ]] || [[ "$TUMOR_PREFIX" =~ ^#.* ]]; then
        continue
    fi

    # 检查是否成功读取了两个样本名称
    if [[ -z "$NORMAL_PREFIX" ]]; then
        echo "警告：发现不完整的样本对行: '$TUMOR_PREFIX'。该行将被跳过。"
        continue
    fi
    
    echo "=================================================="
    echo "Processing pair: Tumor $TUMOR_PREFIX vs Normal $NORMAL_PREFIX"
    
    # 调用 process_tn_pair 函数，传入肿瘤和正常样本前缀
    # 确保 wxs_variant_call.sh 中定义了 process_tn_pair 函数
    process_tn_pair "$TUMOR_PREFIX" "$NORMAL_PREFIX"
    
    echo "Completed processing for $TUMOR_PREFIX"
    echo "=================================================="
done < "$SAMPLE_PAIRS_FILE"

echo "Pipeline finished for all sample pairs listed."
