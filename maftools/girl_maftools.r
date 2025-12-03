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

# boxplot
lapply(sig_genes_a, function(x){
  dt_vaf_gm <- dt_data_300[Hugo_Symbol == x]
  
  pdf(file = paste('./results/maftools/VAF/GILR300/', x, '.pdf', sep = ''), 3, 3)
  pt <- ggboxplot(
    dt_vaf_gm, x = 'relapse', y = "vaf", color = 'relapse', fill = 'relapse', 
    palette = color_phe[unique(as.character(dt_vaf_gm$relapse))],
    alpha = 0.6, add = "jitter", add.params = list(size = 0.05, jitter = 0.2, alpha = 2),
    outlier.shape = NA) +
    theme_few() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    labs(x = x, y = 'VAF')
  print(pt)
  dev.off()
})

# lollipopplot
lapply(gene_int_gilr300_aml[-60], function(x){
  if(any(unique(dt_maf_300[x, on = 'Hugo_Symbol'][['HGVSp_Short']]) != '')){
    print(x)
    pdf(file = paste('./results/maftools/GILR300/lollipopplot/', x, '.pdf'), 4, 4)
    pt <- lollipopPlot(df_maf_300, gene = x, AACol = "HGVSp_Short")
    print(pt)
    dev.off()
  }
})

# 

### GILR 3000 ------
# data input ---
dt_maf_3000 <- read.maf('./data/GILR/GILR_3000.maf')
dt_data_3000 <- dt_maf_3000@data
dt_girl_from_company <- fread('/Users/stead/JS_data/level3_data/wang-lab/M13/data/wes/300_3000_GIL.txt')
intersect(dt_girl_from_company[Tumor_Sample_Barcode == '3000_GIL'][['Hugo_Symbol']], dt_data_3000$Hugo_Symbol)

### DMGs ------
e_vs_l <- mafCompare(m1 = df_maf_b, m2 = df_maf_a, m1Name = 'before_relapse', m2Name = 'after_relapse', minMut = 1)
dt_gene_mut <- e_vs_l$results %>% data.table()

###
dt_maf_b <- fread('./data/J002_raw/Somatic_mutation.maf')
dt_maf_a <- fread('./data/J001_relapse/Somatic_mutation.maf')
dt_maf <- rbind(dt_maf_b, dt_maf_a)
fwrite(dt_maf, file = './data/menin_combine.maf', sep = '\t')
dt_phe <- cbind(c(unique(dt_maf_b$Tumor_Sample_Barcode), unique(dt_maf_a$Tumor_Sample_Barcode)),
                c(rep('before', 3), rep('after', 3))) %>% data.table()
colnames(dt_phe) <- c('Tumor_Sample_Barcode', 'relapse')

dt_maf <- read.maf('./data/menin_combine.maf', clinicalData = dt_phe)
dt_maf@clinical.data[, relapse := factor(relapse, levels = c('before', 'after'))]

### oncoplot ---
# after occur
sig_genes_a <- dt_gene_mut[before_relapse == 0 & after_relapse > 1][['Hugo_Symbol']]
pdf(file = './results/dmgs/menin_top_dmgs_after.pdf', 6, 24)
oncoplot(maf = dt_maf, colors = vc_cols, draw_titv = TRUE, genes = sig_genes_a,
         clinicalFeatures = c('relapse'), removeNonMutated = FALSE,
         showTumorSampleBarcodes = TRUE,
         annotationColor = list(
           relapse = color_phe[names(color_phe) %in% c('before', 'after')],
           sortByAnnotation = FALSE, anno_height = 0.5)
)
dev.off()

# before occur
sig_genes_b <- dt_gene_mut[before_relapse > 1 & after_relapse  == 0][['Hugo_Symbol']]
pdf(file = './results/dmgs/menin_top_dmgs_before.pdf', 6, 24)
oncoplot(maf = dt_maf, colors = vc_cols, draw_titv = TRUE, genes = sig_genes_b,
         clinicalFeatures = c('relapse'), removeNonMutated = FALSE,
         showTumorSampleBarcodes = TRUE,
         annotationColor = list(
           relapse = color_phe[names(color_phe) %in% c('before', 'after')],
           sortByAnnotation = FALSE, anno_height = 0.5)
)
dev.off()

### VAF after ---
dt_data <- dt_maf@data
lapply(sig_genes_a, function(x){
  dt_vaf_gm <- dt_data[Hugo_Symbol == x][dt_phe, on = 'Tumor_Sample_Barcode']
  
  pdf(file = paste('./results/vaf/after/', x, '.pdf', sep = ''), 3, 3)
  pt <- ggboxplot(
    dt_vaf_gm, x = 'relapse', y = "t_AF", color = 'relapse', fill = 'relapse', 
    palette = color_phe[unique(as.character(dt_vaf_gm$relapse))],
    alpha = 0.6, add = "jitter", add.params = list(size = 0.05, jitter = 0.2, alpha = 2),
    outlier.shape = NA) +
    theme_few() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    labs(x = x, y = 'VAF')
  print(pt)
  dev.off()
})

### VAF before ---
lapply(sig_genes_b, function(x){
  dt_vaf_gm <- dt_data[Hugo_Symbol == x][dt_phe, on = 'Tumor_Sample_Barcode']
  
  pdf(file = paste('./results/vaf/before/', x, '.pdf', sep = ''), 3, 3)
  pt <- ggboxplot(
    dt_vaf_gm, x = 'relapse', y = "t_AF", color = 'relapse', fill = 'relapse', 
    palette = color_phe[unique(as.character(dt_vaf_gm$relapse))],
    alpha = 0.6, add = "jitter", add.params = list(size = 0.05, jitter = 0.2, alpha = 2),
    outlier.shape = NA) +
    theme_few() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    labs(x = x, y = 'VAF')
  print(pt)
  dev.off()
})


### lollipopPlot after ---
lapply(sig_genes_a[-c(12, 19, 22, 27, 30, 34, 37, 51, 57, 72)], function(x){
  print(x)
  pdf(file = paste('./results/lollipopplot/after/', x, '.pdf'), 4, 4)
  pt <- lollipopPlot(dt_maf, gene = x, AACol = "AAchange")
  print(pt)
  dev.off()
})

### lollipopPlot before ---
lapply(sig_genes_b[-c(19, 27, 29, 30, 34, 46, 55, 64, 70, 75)], function(x){
  print(which(sig_genes_b == x))
  pdf(file = paste('./results/lollipopplot/before/', x, '.pdf'), 4, 4)
  pt <- lollipopPlot(dt_maf, gene = x, AACol = "AAchange")
  print(pt)
  dev.off()
})

dt_en1 <- dt_data[Hugo_Symbol == 'EN1']
dt_cops6 <- dt_data[Hugo_Symbol == 'COPS6']
dt_rdm1 <- dt_data[Hugo_Symbol == 'RDM1']
dt_men1 <- dt_data[Hugo_Symbol == 'MEN1']

dt_data <- dt_data[Hugo_Symbol %in% c('COPS6', 'MEN1', 'RDM1', 'NECTIN2')]

# relase sample over two samples ---
dt_data <- dt_data[Hugo_Symbol %in% sig_genes_a]
dt_data_f <- dt_data[, c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 
                         'Tumor_Seq_Allele2', 't_AF', 'AAchange'), with = FALSE]
fwrite(dt_data_f, file = './data/relapse_gene_over_2sample_250929.csv')

### relapse and high vaf ---
dt_data_relapse <- dt_data[Hugo_Symbol %in% sig_genes_a]
dt_data_relapse_0.2 <- dt_data_relapse[t_AF > 0.2][, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 
                                                       'Tumor_Seq_Allele2', 't_AF', 'AAchange'), with = FALSE]
fwrite(dt_data_relapse_0.2, file = './data/relapse_gene_vaf_0.2.csv')

# no filter
dt_data_relapse <- dt_data_relapse[, c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 
                                       'Tumor_Seq_Allele2', 't_AF', 'AAchange'), with = FALSE]
fwrite(dt_data_relapse, file = './data/relapse_gene_vaf_all.csv')

### VAF degs ---
dt_mut <- dt_maf@data
dt_mut[Tumor_Sample_Barcode %in% unique(dt_maf_b$Tumor_Sample_Barcode), group := 'before']
dt_mut[Tumor_Sample_Barcode %in% unique(dt_maf_a$Tumor_Sample_Barcode), group := 'after']

# paired vaf up ---
dt_mut[, patient_id := unlist(lapply(strsplit(as.character(Tumor_Sample_Barcode), '_'), '[[', 1))]
dt_mut[, pos_id := paste(patient_id, Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = '_')]
dt_mut_u_before <- dt_mut[group  == 'before']
before_vaf_0 <- setdiff(unique(dt_mut$pos_id), unique(dt_mut_u_before$pos_id))
dt_mut_before_0 <- cbind(unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 1)), 
                         unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 2)), 
                         paste(unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 1)), 
                               unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 2)), 
                               unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 3)),
                               unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 4)),
                               unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 5)),
                               unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 6)), 
                               unlist(lapply(strsplit(before_vaf_0, '_'), '[[', 7)), sep = '_'),
                         'before', 0) %>% data.table()
colnames(dt_mut_before_0) <- c('patient_id', 'Hugo_Symbol', 'pos_id', 'group', 't_AF')
dt_mut_u_before_m <- rbind(dt_mut_u_before[, c('patient_id', 'Hugo_Symbol', 'pos_id', 'group', 't_AF'), with = FALSE],
                           dt_mut_before_0)

# after
dt_mut_u_after <- dt_mut[group  == 'after']
after_vaf_0 <- setdiff(unique(dt_mut$pos_id), unique(dt_mut_u_after$pos_id))
dt_mut_after0 <- cbind(unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 1)), 
                       unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 2)), 
                       paste(unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 1)), 
                             unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 2)), 
                             unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 3)),
                             unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 4)),
                             unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 5)),
                             unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 6)),
                             unlist(lapply(strsplit(after_vaf_0, '_'), '[[', 7)), sep = '_'),
                       'after', 0) %>% data.table()
colnames(dt_mut_after0) <- c('patient_id', 'Hugo_Symbol', 'pos_id', 'group', 't_AF')
dt_mut_u_after_m <- rbind(dt_mut_u_after[, c('patient_id', 'Hugo_Symbol', 'pos_id', 'group', 't_AF'), with = FALSE],
                          dt_mut_after0)

dt_mut_m <- rbind(dt_mut_u_before_m, dt_mut_u_after_m)
dt_mut_mw <- dt_mut_u_before_m[dt_mut_u_after_m, on = 'pos_id']
dt_vaf_genes <- lapply(unique(dt_mut_mw$Hugo_Symbol), function(x){
  dt_mut_mwd <- dt_mut_mw[x, on = 'Hugo_Symbol']
  dt_mut_mwd[, vaf_f := as.numeric(i.t_AF) - as.numeric(t_AF)]
  dt_mut_mwd[, vaf_p := (as.numeric(i.t_AF) - as.numeric(t_AF))/(as.numeric(t_AF) + 0.0001)]
  if(dim(dt_mut_mwd)[1] > 2 & all(dt_mut_mwd$vaf_f > 0) & all(dt_mut_mwd$i.t_AF > 0.01) & all(dt_mut_mwd$vaf_p > 0.1)){
    return(list(
      gene_name = x
    ))
  }
}) %>% rbindlist() %>% data.table()

lapply(dt_vaf_genes$gene_name, function(x){
  print(x)
  dt_mut_u <- dt_mut_m[x, on = 'Hugo_Symbol']
  if(length(unique(dt_mut_u$patient_id)) > 0){
    
    # test
    #dt_pair <- data.frame(before = dt_mut_u[group == 'before'], after = dt_mut_u[group == 'after'])
    pdf(file = paste('./results/vaf/pair/', x, '.pdf', sep = ''), 3, 3)
    pt <- ggpaired(dt_mut_u, x = "group", y = "t_AF", 
                   color = "group", palette = color_phe[unique(dt_mut_u$group)], id = 'pos_id',
                   line.color = "gray", line.size = 0.5, short.panel.labs = FALSE) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
      labs(x = NULL, y = 'VAF', title = x) +
      stat_compare_means(paired = TRUE, hide.ns = TRUE, label = 'p.format')
    print(pt)
    dev.off()
  }
})

# facet
up_genes <- dt_vaf_genes$gene_name
dt_mut_m_up <- dt_mut_m[up_genes, on = 'Hugo_Symbol']
pdf(file = './results/vaf/pair/vaf_up_genes.pdf', 6, 6)
pt <- ggpaired(dt_mut_m_up, x = "group", y = "t_AF", 
               color = "group", palette = color_phe[unique(dt_mut_um$group)], id = 'pos_id',
               line.color = "gray", line.size = 0.5, short.panel.labs = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  facet_wrap(. ~ Hugo_Symbol, scales = 'free') + 
  labs(x = 'group', y = 'VAF')
print(pt)
dev.off()

# heatmap
dt_mut_m_up <- dt_mut_m[dt_vaf_genes$gene_name, on = 'Hugo_Symbol']
dt_mut_m_up[, pos_id2 := paste(unlist(lapply(strsplit(pos_id, '_'), '[[', 2)), 
                               unlist(lapply(strsplit(pos_id, '_'), '[[', 3)),
                               unlist(lapply(strsplit(pos_id, '_'), '[[', 4)),
                               unlist(lapply(strsplit(pos_id, '_'), '[[', 5)),
                               unlist(lapply(strsplit(pos_id, '_'), '[[', 6)),
                               unlist(lapply(strsplit(pos_id, '_'), '[[', 7)), sep = '_')]
dt_mut_m_upw <- dcast(dt_mut_m_up, pos_id2 ~ patient_id + group, value.var = 't_AF')

# clinical annotation
cli_var <- c('group')
dt_phe <- cbind(colnames(dt_mut_m_upw[, !'pos_id2']), c('after', 'before', 'after', 'before', 'after', 'before')) %>% data.table()
colnames(dt_phe) <- c('sample_id', 'group')

pt_bottom <- HeatmapAnnotation(
  df = dt_phe, 
  na_col = 'lightgrey',
  simple_anno_size = unit(0.35, "cm")
)

mat_vaf <- as.matrix(apply(dt_mut_m_upw[, !'pos_id2'], 2, function(x){as.numeric(x)}))
mat_vaf[is.na(mat_vaf)] <- 0
ht = Heatmap(mat_vaf, name = "VAF", col = colorRamp2(c(0, 0.5), c('white', "brown")),
             top_annotation = pt_bottom, column_split = dt_phe[['group']],
             show_row_names = FALSE, show_column_names = TRUE, 
             cluster_columns = FALSE, column_dend_reorder = FALSE, clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2', 
             cluster_rows = TRUE, row_dend_reorder = TRUE, column_names_gp = gpar(fontsize = 16, fontface = "bold"),
             show_row_dend = TRUE, show_column_dend = TRUE, show_heatmap_legend = TRUE) 

pdf(file = './results/vaf/vaf_up_genes.pdf', width = 4, height = 12)
print(ht)
dev.off()

### pathway ---
# Oncogenic Signaling Pathways
pdf(file = './results/pathway/oncogenic_pathway.pdf', 4, 6)
dt_pathways <- OncogenicPathways(maf = dt_maf, pathways = NULL)
dev.off()

# oncoplot
lapply(dt_pathways$Pathway, function(x){
  p_height <- log2(dt_pathways[x, on = 'Pathway'][['n_affected_genes']])
  
  if(x %in% c('RTK-RAS', 'WNT', 'Hippo', 'NOTCH', 'PI3K')){
    pw = 3
    ph = 8
  } else {
    pw = 3
    ph = 3
  }
  pdf(file = paste('./results/pathway/', x, '.pdf', sep = ''), pw, ph)
  PlotOncogenicPathways(maf = dt_maf, pathways = x, showTumorSampleBarcodes = TRUE, 
                        removeNonMutated = FALSE, fullPathway = TRUE,
                        sampleOrder = c('ZSJ_screen', 'MLJ_screen', 'XXP_screen', 
                                        'ZSJ_20241113', 'MLJ_20250124', 'XXP_20250103'))
  dev.off()
})

