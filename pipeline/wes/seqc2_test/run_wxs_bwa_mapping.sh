#!/bin/bash

# ------------------------------------------------
# 全局配置变量 (Global Configuration Variables)
# ------------------------------------------------
data_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq
results_base=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/results
log_dir=/mnt/pgx_src_data_pool_4/seqc2/wes/fastq/gatk_test/log
fasta=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/gatk_resource/Homo_sapiens_assembly38.fasta

# Input file suffixes
R1_suffix="_1.fq.gz"
R2_suffix="_2.fq.gz"

# Execution environments (Singularity images)
IMG_BWA=/mnt/pgx_src_data_pool_4/reference/images/dnaseq-sv-workflow-v0.1.img
IMG_SAMTOOLS=/mnt/pgx_src_data_pool_4/reference/images/rnaseq-exp-workflow-v0.1.img
IMG_GATK=/mnt/pgx_src_data_pool_4/reference/images/gatk-v4.6.2.0.img

# Number of threads (CPU cores)
nt=16

# Sequencing platform
pl="ILLUMINA" 

# 确保日志目录存在
mkdir -p "$log_dir"

# ------------------------------------------------
# 引入核心处理函数 (Source the core processing function)
# ------------------------------------------------
source ./wxs_bwa_mapping.sh

# ------------------------------------------------
# 定义需要处理的样本列表 (Define the list of samples to process)
# ------------------------------------------------
SAMPLES_TO_PROCESS=(
    "WES_IL_N_1"
    "WES_IL_T_1" 
    # "在此处添加更多样本名称，例如 Sample_002"
    # "ThirdSample_ABC_L003"
)

# ------------------------------------------------
# 循环处理所有样本 (Loop through all samples)
# ------------------------------------------------
for sample_name in "${SAMPLES_TO_PROCESS[@]}"; do
    echo "=================================================="
    echo "Initiating analysis for: $sample_name"
    
    # 调用 process_sample 函数，并将当前样本名称作为参数传入
    process_sample "$sample_name"
    
    echo "Completed analysis for: $sample_name"
    echo "=================================================="
done

echo "Pipeline finished for all samples listed."
