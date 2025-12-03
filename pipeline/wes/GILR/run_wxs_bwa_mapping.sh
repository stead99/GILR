#!/bin/bash

# ------------------------------------------------
# 全局配置变量 (Global Configuration Variables)
# ------------------------------------------------
data_dir=/mnt/pgx_src_data_pool_4/GILR/data/WES
results_base=/mnt/pgx_src_data_pool_4/GILR/results/WES
log_dir=/mnt/pgx_src_data_pool_4/GILR/code/wes_variant_call/log
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
# 从外部文件读取样本列表 (Read sample list from an external file)
# ------------------------------------------------

# 定义样本列表文件的路径和名称
# 默认使用当前目录下的 samples_list.txt 文件
SAMPLE_LIST_FILE="./sample_girl3000.txt"

# 检查样本列表文件是否存在
if [[ ! -f "$SAMPLE_LIST_FILE" ]]; then
    echo "错误：找不到样本列表文件 '$SAMPLE_LIST_FILE'。"
    echo "请创建一个名为 samples_list.txt 的文件，并在其中每行放置一个样本名称。"
    exit 1
fi

# 使用 mapfile 或 while read 循环读取文件内容
# mapfile -t SAMPLES_TO_PROCESS < "$SAMPLE_LIST_FILE"
# 使用 while read 是更具兼容性的方式
echo "正在从 $SAMPLE_LIST_FILE 读取样本列表..."

# ------------------------------------------------
# 循环处理所有样本 (Loop through all samples)
# ------------------------------------------------
while IFS= read -r sample_name; do
    # 跳过空行或注释行（以 # 开头）
    if [[ -z "$sample_name" ]] || [[ "$sample_name" =~ ^#.* ]]; then
        continue
    fi
    
    echo "=================================================="
    echo "Initiating analysis for: $sample_name"
    
    # 调用 process_sample 函数，并将当前样本名称作为参数传入
    # 确保 wxs_bwa_mapping.sh 中定义了 process_sample 函数
    process_sample "$sample_name"
    
    echo "Completed analysis for: $sample_name"
    echo "=================================================="
done < "$SAMPLE_LIST_FILE"

echo "Pipeline finished for all samples listed."
