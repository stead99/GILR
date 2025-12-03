library(data.table)
library(maftools)
source('./code/config_files.R')

### 300 ---
file_names <- c('CTRL', 'GIL')
dt_300_gs <- do.call(rbind, lapply(file_names, function(x){
  # snv 
  dt_snv_u <- fread(paste('./data/wes/M13_300/', 'M13_', x, '.samtools.snp.annovar.hg19_multianno.xls.gz', sep = ''))
  dt_snv_u_coding <- dt_snv_u[c('exonic'), on = 'Func']
  dt_snv_u_coding[, Variant_Type := 'SNP']
  
  # indel
  dt_indel_u <- fread(paste('./data/wes/M13_300/', 'M13_', x, '.samtools.indel.annovar.hg19_multianno.xls.gz', sep = ''))
  dt_indel_u_coding <- dt_indel_u[c('exonic'), on = 'Func']
  dt_indel_u_coding[ExonicFunc %in% c('frameshift deletion', 'nonframeshift deletion'), Variant_Type := 'DEL']
  dt_indel_u_coding[ExonicFunc %in% c('frameshift insertion', 'nonframeshift insertion'), Variant_Type := 'INS']
  
  dt_var <- rbind(dt_snv_u_coding, dt_indel_u_coding)
  colnames(dt_var)[57] <- 'type'
  dt_var[, Tumor_Sample_Barcode := paste('300_', x, sep = '')]
  dt_var[, var_id := paste(CHROM, POS, REF, ALT, sep = '_')]
  return(dt_var)
  })
  )

dt_300 <- dt_300_gs[Tumor_Sample_Barcode == '300_GIL'][!unique(dt_300_gs[Tumor_Sample_Barcode == '300_CTRL'][['var_id']]), on = 'var_id']
dt_300[, Center := 'Novogene'][, NCBI_Build := 'GRCh37'][, Start_position := POS][, End_Position := POS]
dt_300[, Strand := ''][, Tumor_Seq_Allele1 := REF]
dt_300_d <- dt_300[, c('GeneName', 'Gene', 'Center', 'NCBI_Build', 'CHROM', 'Start_position', 'End_Position', 'Strand', 
                       'ExonicFunc', 'Variant_Type', 'REF', 'Tumor_Seq_Allele1', 'ALT', 'Tumor_Sample_Barcode', 'AAChange'), with = FALSE]
colnames(dt_300_d) <- c('Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_position', 'End_Position',
                        'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                        'Tumor_Sample_Barcode', 'AAChange')
dt_300_dd <- dt_300_d[!c('synonymous SNV', 'unknown', 'stopgain'), on = 'Variant_Classification']
dt_300_dd[, Chromosome := paste('chr', Chromosome)]
dt_300_dd[Variant_Classification == 'missense SNV', Variant_Classification := 'Missense_Mutation']
dt_300_dd[Variant_Classification == 'frameshift deletion', Variant_Classification := 'Frame_Shift_Del']
dt_300_dd[Variant_Classification == 'frameshift insertion', Variant_Classification := 'Frame_Shift_Ins']
dt_300_dd[Variant_Classification == 'nonframeshift insertion', Variant_Classification := 'In_Frame_Ins']
dt_300_dd[Variant_Classification == 'nonframeshift deletion', Variant_Classification := 'In_Frame_Del']

### 3000
file_names <- c('gilr', 'WT')
dt_3000 <- lapply(file_names, function(x){
  # snv 
  dt_snv_u <- fread(paste('./data/wes/M13_3000/', x, '.bcftools.snp.annovar.hg38_multianno.xls.gz', sep = ''))
  dt_snv_u_coding <- dt_snv_u[c('exonic'), on = 'Func']
  dt_snv_u_coding[, Tumor_]
  
  # indel
  dt_indel_u <- fread(paste('./data/wes/M13_3000/', x, '.bcftools.indel.annovar.hg38_multianno.xls.gz', sep = ''))
  dt_indel_u_coding <- dt_indel_u[c('exonic'), on = 'Func']
})

###
dt_3000 <- fread('./data/wes/M13_3000/Somatic_mutation.maf')
dt_3000_exon <- dt_3000[c('Missense_Mutation'), on = 'Variant_Classification']
dt_3000_exon[, Tumor_Sample_Barcode := '3000_GIL']

col_int <- intersect(colnames(dt_300_dd), colnames(dt_3000_exon))
dt_m <- rbind(dt_3000_exon[, col_int, with = FALSE],
              dt_300_dd[, col_int, with = FALSE])
fwrite(dt_m, file = './data/wes/300_3000_GIL.txt', sep = '\t')
fwrite(dt_m, file = './data/wes/300_3000_GIL.maf', sep = '\t')

### maftools ---
df_maf <- read.maf('./data/wes/300_3000_GIL.maf')
dt_phe <- cbind(unique(df_maf@data$Tumor_Sample_Barcode), unique(df_maf@data$Tumor_Sample_Barcode)) %>% data.table()
colnames(dt_phe) <- c('Tumor_Sample_Barcode', 'group')
df_maf <- read.maf('./data/wes/300_3000_GIL.maf', clinicalData = dt_phe)

pdf(file = './results/300_3000_GIL.pdf', 6, 24)
oncoplot(maf = df_maf, colors = vc_cols, draw_titv = TRUE,
         removeNonMutated = FALSE, top = 50,
         showTumorSampleBarcodes = TRUE,
         annotationColor = list(
           group = color_phe[names(color_phe) %in% c('300_GIL', '3000_GIL')],
           sortByAnnotation = FALSE, anno_height = 0.5)
         )
dev.off()
