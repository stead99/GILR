library(data.table)
library(ggplot2)
library(VennDiagram)

### SEQC2 high confidence region ------
seqc2_dir <- '/Users/stead/JS_data/source_tree/project_code/public_github/somaticseq-seqc2_v1.2/ftp_data/release/v1.2/'

v_type = 'snv'

# seqc2 high snv
if(v_type == 'snv'){
  hc_path <- paste(seqc2_dir, 'high-confidence_sSNV_in_HC_regions_v1.2.vcf', sep = '')
} else if (v_type == 'indel'){
  hc_path <- paste(seqc2_dir, 'high-confidence_sINDEL_in_HC_regions_v1.2.vcf', sep = '')
}

snv_pos_h <- grep('#CHROM', readLines(hc_path), value = FALSE)
dt_snv_h <- fread(hc_path, skip = snv_pos_h - 1, sep = '\t', header = TRUE)
setnames(dt_snv_h, '#CHROM', 'CHR')
dt_snv_h[, pos_num := paste(CHR, POS, REF, ALT, sep = '_')]

### compare with high confidence ------
# filter position
get_pos <- function(data_dir, filename){
  # test file
  vcf_file <- paste(data_dir, filename, sep = '')
  head_pos <- grep('#CHROM', readLines(vcf_file), value = FALSE)
  dt_v <- fread(vcf_file, skip = head_pos - 1, sep = '\t', header = TRUE)
  
  if(dim(dt_v[FILTER  == 'PASS'])[1] > 1){
    dt_v <- dt_v[FILTER  == 'PASS']
  } 
  
  setnames(dt_v, '#CHROM', 'CHR')
  dt_v[, pos_num := paste(CHR, POS, REF, ALT, sep = '_')]
  return(dt_v)
}


# vcf position compare recall, precision and f-score
vcf_compare <- function(dt_hc, dt_test, test_name, res_name){
  veen_hc <- dt_hc$pos_num
  veen_test <- dt_test$pos_num
  
  # veen plot
  veen_list <- list(veen_hc, veen_test)
  names(veen_list) <- c('seqc2_hc', test_name)
  pt_veen <- venn.diagram(veen_list,
                          filename = NULL, fill = c('peru', 'deepskyblue2'), 
                          alpha = c(0.3, 0.3), cex = 2,
                          col = c('peru', 'deepskyblue2'), lwd = 2,
                          cat.col  = c('peru', 'deepskyblue2'), cat.cex = 2,  cat.fontface = "bold",
                          cat.default.pos = "outer")
  pdf(file = paste('./results/seqc2/', res_name, '.pdf', sep = ''), 6, 6)
  grid.draw(pt_veen)
  dev.off()
  
  # recall, precion and f-score
  recall = length(intersect(veen_hc, veen_test))/(length(intersect(veen_hc, veen_test)) + length(setdiff(veen_hc, veen_test)))
  precision = length(intersect(veen_hc, veen_test))/(length(intersect(veen_hc, veen_test)) + length(setdiff(veen_test, veen_hc)))
  fscore = 2*(precision * recall)/(precision + recall)
  total_count = length(veen_test)
  
  return(
    list(total_count = total_count, 
         fscore = fscore, 
         recall  = recall,
         precision = precision)
  )
}

### test with hc regeion ---
dt_test <- get_pos(data_dir = './data/seqc2/',
                   filename = 'WES_IL_T_1.mutect2.PASS.snv.hc.vcf.gz')
dt_performamce <- vcf_compare(dt_hc = dt_snv_h, dt_test = dt_test, 
                              test_name = 'three_consensus_il_1',
                              res_name = 'seqc2_three_consensus_il_1_snv')

### IL1-3 tnseq and tnscope intersection
dt_int_p <- lapply(1, function(x){
  
  # tnseq
  dt_tnseq <- get_pos(data_dir = paste('./data/seqc2/', v_type, '/', sep = ''),
                      filename = paste('WGS_IL_T_', x, '_WGS_IL_N_', x, paste('.TNseq.PASS.', v_type, '.hc.recode.vcf', sep = ''), sep = ''))
  # tnscope
  dt_tnscope <- get_pos(data_dir =  paste('./data/pipeline_test/variant_call/BAM2VCF/', v_type, '/', sep = ''), 
                        filename = paste('WGS_IL_T_', x, '_WGS_IL_N_', x, paste('.TNscope.PASS.', v_type, '.hc.recode.vcf', sep = ''), sep = ''))
  
  # strelka
  dt_strelka <-  get_pos(data_dir =  paste('./data/pipeline_test/variant_call/BAM2VCF/', v_type, '/', sep = ''), 
                         filename = paste('WGS_IL_T_', x, '_WGS_IL_N_', x, paste('.strelka.PASS.', v_type, '.hc.recode.nobqsr.vcf', sep = ''), sep = ''))
  
  pos_num <- names(which(table(c(dt_tnseq$pos_num, dt_tnscope$pos_num, dt_strelka$pos_num)) >= 2))
  int_sofware <- table(c(dt_tnseq$pos_num, dt_tnscope$pos_num, dt_strelka$pos_num))[pos_num]
  
  # intersection 
  dt_int <- cbind(pos_num, int_sofware) %>% data.table()
  
  # compare performance
  performamce_list_tnscope <- vcf_compare(dt_hc = dt_snv_h, dt_test = dt_int, 
                                          test_name = paste('tnseq_scope_strelka2_il_', v_type, x, sep = ''),
                                          res_name = paste('seqc2_tnseq_scope_strelka2_il_', v_type, x, sep = ''))
  performamce_list_tnscope$samples <- paste('WGS_IL_T_', x, sep = '')
  performamce_list_tnscope$caller <- 'TNseq_TNScope_strelka2'
  performamce_list_tnscope$type <- 'Fresh'
  
  return(performamce_list_tnscope)
}) %>% rbindlist() %>% data.table()


### snv performance ------
# performance
dt_p <- fread('./data/pipeline_test/performance/performance_test.txt')
dt_test_po <- dt_int_p[, colnames(dt_p), with = FALSE]
dt_p_a <- rbind(dt_p, dt_test_po)

pt <- ggplot(dt_p_a, aes(x = precision, y = recall, shape = type, color = caller)) + 
  geom_point(size = 2.5)+
  scale_shape_manual(values = c(16, 17))+
  theme_classic()+
  xlab('Precision') +
  ylab('Recall') +
  scale_color_brewer(palette = "Set1") +
  theme_few() +
  theme(axis.text.x = element_text(size = 11, color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11, color='black'))+
  theme(axis.title.y = element_text(size = 11, color='black'))+
  theme(axis.title.x = element_text(size = 11, color='black'))

pdf('./results/pipeline_test/performance_tnseq_scope_il.pdf',height = 3, width = 5.5)
plot(pt)
dev.off()


# barplot
pt_bar <- ggplot(data = dt_p_a[type == 'Fresh'], aes(x = reorder(caller, fscore), y = fscore, fill = samples)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_brewer(palette = "Blues") +
  xlab('Caller') +
  ylab('F-Score') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 


pdf('./results/pipeline_test/performance_tnseq_scope_il_bar.pdf',height = 4, width = 4)
plot(pt_bar)
dev.off()

