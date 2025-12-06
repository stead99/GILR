#!/bin/bash

# ======================================================
# 1. 定义全局配置和路径 (Global Configuration)
# ======================================================
log_dir=/mnt/pgx_src_data_pool_4/GILR/code/wes_variant_call/variant_call
data_dir=/mnt/pgx_src_data_pool_4/GILR/results/WES_bobed/Mutect2
vep_dir=/mnt/pgx_src_data_pool_4/GILR/results/WES_bobed/VEP
anno_dir=/mnt/pgx_src_data_pool_4/GILR/results/WES_bobed/annovar
fasta=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/sentieon/GRCh38.d1.vd1.fa
cache=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/ensembl_vep
annovar_database=/mnt/pgx_src_data_pool_4/reference/human/GRCh38/annovar_database
nt=6
IMG_VEP=/mnt/pgx_src_data_pool_4/reference/images/vep-v104.0.img
IMG_annovar=/mnt/pgx_src_data_pool_4/reference/images/annovar-v201910.img

# 确保必要的目录存在
mkdir -p "${vep_dir}"
mkdir -p "${anno_dir}"
mkdir -p "${log_dir}/log"

# ======================================================
# 2. 指定外部样本列表文件名 (Specify External Samples List File)
# ======================================================
# 将文件名直接硬编码在这里，不再需要命令行参数
SAMPLES_FILE="sample_anno_gilr.txt"

if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file '$SAMPLES_FILE' not found in the current directory."
    exit 1
fi

# ======================================================
# 3. 批量处理循环 (Batch Processing Loop)
# ======================================================
echo "开始批量处理样本 (从文件 ${SAMPLES_FILE} 读取)..."

# 使用 while 循环直接读取外部文件内容
while IFS=$' \t' read -r vcf tumor_id normal_id || [ -n "$vcf" ]; do
    
    # 跳过注释行 (#) 和空行
    if [[ "$vcf" =~ ^#.* ]] || [[ -z "$vcf" ]]; then
        continue
    fi

    echo "--- 正在处理样本: ${tumor_id} (VCF: ${vcf}) ---"

    # 设置当前样本的日志文件
    logfile="$tumor_id.$(date +"%Y%m%d%H%M").log"
    # 将标准错误输出重定向到该日志文件
    exec 2> "${log_dir}/log/$logfile"
    set -x # 在日志中开启命令追踪

    # ******************************************
    # 1. VEP 注释流程
    # ******************************************
    start_time_mod=$(date +%s)

    singularity exec "${IMG_VEP}" \
    perl /opt/vep/ensembl-vep/vep --format vcf --vcf \
        --assembly GRCh38 \
        --species homo_sapiens_merged \
        --everything --af_exac \
        --offline \
        --cache --dir_cache "${cache}" \
        --fasta "${fasta}" \
        --sift b \
        --polyphen b \
        --input_file "${data_dir}/${vcf}" --output_file "${vep_dir}/${tumor_id}.vep.vcf"

    # ******************************************
    # 2. VCF to MAF 转换流程
    # ******************************************
    singularity exec "${IMG_VEP}" \
    perl /opt/mskcc-vcf2maf/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf "${vep_dir}/${tumor_id}.vep.vcf" \
        --output-maf "${vep_dir}/${tumor_id}.maf" \
        --tumor-id "${tumor_id}" \
        --normal-id "${normal_id}" \
        --ref-fasta "${fasta}" \
        --ncbi-build GRCh38 \
        --species homo_sapiens_merged \
        --vep-fork "$nt"

    # ******************************************
    # 1. annovar
    # ******************************************
    start_time_mod=$(date +%s)
    singularity exec  "${IMG_annovar}" \
    table_annovar.pl "${data_dir}/${vcf}" \
        ${annovar_database} -buildver hg38 \
        -out ${anno_dir}/${tumor_id} -remove \
        -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,esp6500siv2_all,ALL.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EAS.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08,avsnp147,dbnsfp33a,clinvar_20210501,gnomad_genome,dbscsnv11,dbnsfp31a_interpro \
        -operation g,g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f \
        -nastring . -vcfinput -thread $nt
        
    set +x # 关闭命令追踪

    end_time_mod=$(date +%s)
    if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" "$start_time_mod"); else start_date=$(date -d "@$start_time_mod"); fi
    if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" "$end_time_mod"); else end_date=$(date -d "@$end_time_mod"); fi
    echo "Module VEP_annovar Started: $start_date; Ended: $end_date; Elapsed time: $(($end_time_mod - $start_time_mod)) sec">> "${log_dir}/log/$logfile"
    
    echo "--- 样本 ${tumor_id} 处理完成 ---"

# 循环从硬编码的文件名变量 $SAMPLES_FILE 读取
done < "$SAMPLES_FILE"

echo "所有样本批量处理完毕。"
