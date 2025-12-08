library(maftools)
library(dplyr)
library(data.table)
source('./code/config_files.R')

dt_aml_genes <- fread('./data/AML/Mutated_Genes_cbio_hg19.txt')
dt_aml_genes[, freq_n := as.numeric(gsub('%', '', Freq))]
high_aml_genes <- dt_aml_genes[freq_n >= 0.5][['Gene']]
dt_bed <- fread('./data/bed/Exome-Agilent_V6.bed')

### GILR 300 ------
# filter no bed region point ---
dt_maf_300 <- fread('./data/GILR/nobed/M13_GIL_FKDO230138540-1A_HH2JTDSX5_L2.maf')
dt_maf_300d <- do.call(rbind, lapply(seq(dim(dt_maf_300)[1]), function(x){
  chr_n = dt_maf_300[x, ][['Chromosome']]
  start_n = dt_maf_300[x, ][['Start_Position']]
  end_n = dt_maf_300[x, ][['End_Position']]
  gene_n = dt_maf_300[x, ][['Hugo_Symbol']]
  dt_bed_d <- dt_bed[V1 == chr_n & V4 == gene_n][V2 < start_n & V3 > end_n]
  if(dim(dt_bed_d)[1] > 0){
    return(dt_maf_300[x, ])
  }
  })
)
dt_maf_300dm <- rbind(dt_maf_300d, dt_maf_300['FLT3', on = 'Hugo_Symbol'][Start_Position == 28027222])
fwrite(dt_maf_300dm, file = './data/GILR/nobed/GILR_300_filterd_with_bed.maf', sep = '\t')

### data input ---
vc_new = c("Frame_Shift_Del", "Frame_Shift_Ins",  
           "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", 
           "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", 
           "Missense_Mutation")
df_maf_300 <- read.maf('./data/GILR/nobed/GILR_300_filterd_with_bed.maf')

# oncoplot 
dt_data_300 <- df_maf_300@data
dt_data_300[, vaf := t_alt_count/t_depth]
gene_int_gilr300_aml <- intersect(high_aml_genes, dt_data_300$Hugo_Symbol)
pdf(file = './results/maftools/GILR300/GIL_300.pdf', 3, 12)
oncoplot(maf = df_maf_300, colors = vc_cols, draw_titv = TRUE,
         removeNonMutated = FALSE, genes = gene_int_gilr300_aml,
         showTumorSampleBarcodes = TRUE,
         annotationColor = list(
           group = color_phe[names(color_phe) %in% c('300_GIL', '3000_GIL')],
           sortByAnnotation = FALSE, anno_height = 0.5)
)
dev.off()

### VAF ---
# lollipopplot
lapply(gene_int_gilr300_aml[-6], function(x){
  if(any(unique(dt_maf_300[x, on = 'Hugo_Symbol'][['HGVSp_Short']]) != '')){
    print(x)
    pdf(file = paste('./results/maftools/GILR300/lollipopplot/', x, '.pdf', sep = ''), 4, 4)
    pt <- lollipopPlot(df_maf_300, gene = x, AACol = "HGVSp_Short")
    print(pt)
    dev.off()
  }
})

# gene vaf 
dt_data_300d <- dt_data_300[gene_int_gilr300_aml, on = 'Hugo_Symbol']
pdf(file = './results/maftools/GILR300/vaf/vaf_do_plot_gilr300.pdf', 3, 3)
pt <- ggdotchart(dt_data_300d, x = "Hugo_Symbol", y = "vaf",
                 color = "Red",                                # Color by groups
                 sorting = "descending",                        # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 rotate = TRUE,     
                 dot.size = 3,
                 ggtheme = theme_pubr())
print(pt)
dev.off()

### output GILR 300 data table
dt_data_girl300 <- dt_data_300[, colnames(dt_data_300)[c(1:2, 4:7, 9:18, 35:45, 115)], with = FALSE]
fwrite(dt_data_girl300, file = './results/maftools/GILR300/data/gilr300_somatic_mutation.txt', sep ='\t')

### GILR 3000 ------
# filter no bed region point ---
dt_maf_3000 <- fread('./data/GILR/nobed/GILR_3000.maf')
dt_maf_3000d <- do.call(rbind, lapply(seq(dim(dt_maf_3000)[1]), function(x){
  chr_n = dt_maf_3000[x, ][['Chromosome']]
  start_n = dt_maf_3000[x, ][['Start_Position']]
  end_n = dt_maf_3000[x, ][['End_Position']]
  gene_n = dt_maf_3000[x, ][['Hugo_Symbol']]
  dt_bed_d <- dt_bed[V1 == chr_n & V4 == gene_n][V2 < start_n & V3 > end_n]
  if(dim(dt_bed_d)[1] > 0){
    return(dt_maf_3000[x, ])
  }
})
)
dt_maf_3000dm <- rbind(dt_maf_3000d, dt_maf_3000['FLT3', on = 'Hugo_Symbol'][Start_Position == 28027222])
fwrite(dt_maf_3000dm, file = './data/GILR/nobed/GILR_3000_filterd_with_bed.maf', sep = '\t')

# data input ---
df_maf_3000 <- read.maf('./data/GILR/nobed/GILR_3000_filterd_with_bed.maf')

# compared with company ---
dt_data_3000 <- df_maf_3000@data
dt_girl_from_company <- fread('/Users/stead/JS_data/level3_data/wang-lab/M13/data/wes/300_3000_GIL.txt')
intersect(unique(dt_girl_from_company[Tumor_Sample_Barcode == '3000_GIL'][['Hugo_Symbol']]), unique(dt_data_3000$Hugo_Symbol))

# oncoplot
dt_data_3000[, vaf := t_alt_count/t_depth]
gene_int_gilr3000_aml <- intersect(high_aml_genes, dt_data_3000$Hugo_Symbol)
pdf(file = './results/maftools/GILR3000/GIL_3000.pdf', 3, 12)
oncoplot(maf = df_maf_3000, colors = vc_cols, draw_titv = TRUE,
         removeNonMutated = FALSE, genes = gene_int_gilr3000_aml,
         showTumorSampleBarcodes = TRUE,
         annotationColor = list(
           group = color_phe[names(color_phe) %in% c('300_GIL', '3000_GIL')],
           sortByAnnotation = FALSE, anno_height = 0.5)
)
dev.off()

### VAF ---
# lollipopplot
lapply(gene_int_gilr3000_aml[-c(12, 24)], function(x){
  print(x)
  if(any(unique(dt_data_3000[x, on = 'Hugo_Symbol'][['HGVSp_Short']]) != '')){
    print(x)
    pdf(file = paste('./results/maftools/GILR3000/lollipopplot/', x, '.pdf', sep = ''), 4, 4)
    pt <- lollipopPlot(df_maf_3000, gene = x, AACol = "HGVSp_Short")
    print(pt)
    dev.off()
  }
})

# gene vaf 
dt_data_3000d <- dt_data_3000[gene_int_gilr3000_aml, on = 'Hugo_Symbol']
pdf(file = './results/maftools/GILR3000/vaf_do_plot_3000.pdf', 3, 6)
pt <- ggdotchart(dt_data_3000d, x = "Hugo_Symbol", y = "vaf",
                 color = "Red",                                # Color by groups
                 sorting = "descending",                        # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 rotate = TRUE,     
                 dot.size = 3,
                 ggtheme = theme_pubr())
print(pt)
dev.off()
