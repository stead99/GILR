library(maftools)
library(dplyr)
library(data.table)
source('./code/config_files.R')

dt_aml_genes <- fread('./data/AML/Target_GDC_Mutated_Genes.txt')

### GILR 300 ------
dt_maf_300 <- fread('./data/GILR/M13_GIL_FKDO230138540-1A_HH2JTDSX5_L2.maf')

### data input ---
vc_new = c("Frame_Shift_Del", "Frame_Shift_Ins",  
           "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", 
           "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", 
           "Missense_Mutation", 'Intron')
df_maf_300 <- read.maf('./data/GILR/M13_GIL_FKDO230138540-1A_HH2JTDSX5_L2.maf', vc_nonSyn = vc_new)
dt_data_300 <- df_maf_300@data
gene_int_gilr300_aml <- intersect(dt_aml_genes$Gene, dt_data_300$Hugo_Symbol)
pdf(file = './results/maftools/GIL_300.pdf', 3, 12)
oncoplot(maf = df_maf_300, colors = vc_cols, draw_titv = TRUE,
         removeNonMutated = FALSE, genes = c('FLT3', 'STAG2', 'CBL', 'FLG', 'USH2A', 'ATM', 'ARID1A', 'ARID1B', 'APC', 'KMT2C'),
         showTumorSampleBarcodes = TRUE,
         annotationColor = list(
           group = color_phe[names(color_phe) %in% c('300_GIL', '3000_GIL')],
           sortByAnnotation = FALSE, anno_height = 0.5)
)
dev.off()

### VAF ---
# distribution
dt_data_300[, vaf := t_alt_count/t_depth]

# lollipopplot
lapply(gene_int_gilr300_aml[-60], function(x){
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
pdf(file = './results/vaf_do_plot.pdf', 3, 10)
pt <- ggdotchart(dt_data_300d, x = "Hugo_Symbol", y = "vaf",
                 color = "Red",                                # Color by groups
                 sorting = "descending",                        # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 rotate = TRUE,     
                 dot.size = 3,
                 ggtheme = theme_pubr())
print(pt)
dev.off()

### GILR 3000 ------
# data input ---
df_maf_3000 <- read.maf('./data/GILR/GILR_3000.maf', , vc_nonSyn = vc_new)
dt_data_3000 <- df_maf_3000@data
dt_girl_from_company <- fread('/Users/stead/JS_data/level3_data/wang-lab/M13/data/wes/300_3000_GIL.txt')
intersect(unique(dt_girl_from_company[Tumor_Sample_Barcode == '3000_GIL'][['Hugo_Symbol']]), unique(dt_data_3000$Hugo_Symbol))

dt_data_3000 <- df_maf_3000@data
gene_int_gilr3000_aml <- intersect(dt_aml_genes$Gene, dt_data_3000$Hugo_Symbol)
pdf(file = './results/maftools/GIL_3000.pdf', 3, 12)
oncoplot(maf = df_maf_3000, colors = vc_cols, draw_titv = TRUE,
         removeNonMutated = FALSE, genes = c('FLT3', 'STAG2', 'CBL', 'FLG', 'USH2A', 'ATM', 'ARID1A', 'ARID1B', 'APC', 'KMT2C'),
         showTumorSampleBarcodes = TRUE,
         annotationColor = list(
           group = color_phe[names(color_phe) %in% c('300_GIL', '3000_GIL')],
           sortByAnnotation = FALSE, anno_height = 0.5)
)
dev.off()

### VAF ---
# distribution
dt_data_3000[, vaf := t_alt_count/t_depth]

# lollipopplot
lapply(gene_int_gilr3000_aml, function(x){
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
pdf(file = './results/vaf_do_plot_3000.pdf', 3, 15)
pt <- ggdotchart(dt_data_3000d, x = "Hugo_Symbol", y = "vaf",
                 color = "Red",                                # Color by groups
                 sorting = "descending",                        # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 rotate = TRUE,     
                 dot.size = 3,
                 ggtheme = theme_pubr())
print(pt)
dev.off()
