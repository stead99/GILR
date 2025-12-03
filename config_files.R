library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)
library(ggpubr)
library(rstatix)
library(ggthemes)

### mutation type colors ------
vc_cols = c('#1F78B4', '#FB9A99', '#B2DF8A', '#A6CEE3', '#33A02C', '#FF7F00', 
            '#FDBF6F', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', '#E31A1C',  
            '#9ECAE1', "firebrick3", "deepskyblue2", 'plum3', 'plum3', '#FD8D3C', 
            '#1F78B4', '#33A02C', '#6A3D9A', '#E31A1C','white', 
            '#1F78B4', '#A6CEE3', '#B2DF8A', '#33A02C', 'deeppink4', '#FF7F00', '#6A3D9A', '#1F78B4', '#1F78B4', '#1F78B4', 'gray33', # NCCN-hotpot
            '#33A02C', '#1F78B4', '#FF7F00', '#E31A1C', # ecDNA
            "deeppink4", "orange3", "purple3", "gray33", "forestgreen") # sv
names(vc_cols) = c(
  'Missense_Mutation',
  'Frame_Shift_Ins',
  'Nonsense_Mutation',
  'Frame_Shift_Del',
  'Multi_Hit',
  'In_Frame_Del',
  'In_Frame_Ins',
  'Splice_Site',
  'Translation_Start_Site',
  'Truncating',
  'Inframe_Mutation',
  'Promoter_Mutation',
  'Nonstop_Mutation',
  'Amp',
  'Del',
  'Yes',
  'Fusion',
  'pos',
  'clone',
  'early_clone',
  'late_clone',
  'subclone',
   NA, 
  'L858R', 
  'L861Q',
  'S768I', 
  'G719X', 
  'Exon19_Del', 
  'Exon20_Ins', 
  'T790M', 
  'G12C', 
  'V600E',
  'Exon20', 
  'no_hot', 
  'Linear',
  'Complex',
  'BFB',
  'ecDNA',
  'DEL',
  'DUP',
  'INV',
  'CTX', 
  'INS')

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#F8F8F8", col = NA))
  },
  # NA
  miss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "white", col = NA))
  } ,
  # Missense_Mutation
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Missense_Mutation"], col = NA))
  },
  # Nonsense_Mutation
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Nonsense_Mutation"], col = NA))
  },
  # Splice_Site
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Splice_Site"], col = NA))
  },
  # Frame_Shift_Del
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Frame_Shift_Del"], col = NA))
  },
  # Frame_Shift_Ins
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Frame_Shift_Ins"], col = NA))
  },
  # In_Frame_Del
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["In_Frame_Del"], col = NA))
  },
  # In_Frame_Ins
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["In_Frame_Ins"], col = NA))
  },
  # Multi_Hit
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Multi_Hit"], col = NA))
  }, 
  # Translation_Start_Site
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Translation_Start_Site"], col = NA))
  },
  # Truncating
  Truncating = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Truncating"], col = NA))
  },
  # Inframe_Mutation
  Inframe_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Inframe_Mutation"], col = NA))
  },
  # Promoter_Mutation
  Promoter_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Promoter_Mutation"], col = NA))
  },
  # Nonstop_Mutation
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Nonstop_Mutation"], col = NA))
  },
  # Amplification
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Amp"], col = NA))
    },
  # Deletion
  Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = vc_cols["Del"], col = NA))
    },
  # Gene Fusion Yes
  Yes = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Yes"], col = NA))
  },
  # Gene Fusion Fusion
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["Fusion"], col = NA))
  },
  # NCCN target mutation pos or neg
  pos = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["pos"], col = NA))
  },
  # clone
  clone = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["clone"], col = NA))
  },
  # early_clone
  early_clone = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["early_clone"], col = NA))
  },
  # late_clone
  late_clone = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = vc_cols["late_clone"], col = NA))
  },
  # subclone
  subclone = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = vc_cols["subclone"], col = NA))
  },
  # L858R
  L858R = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["L858R"], col = NA))
    },
  # L861Q
    L861Q = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["L861Q"], col = NA))
  },
  # S768I
  S768I = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["S768I"], col = NA))
  },
  # G719X
  G719X = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["G719X"], col = NA))
  },
  # Exon19_Del
  Exon19_Del = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["Exon19_Del"], col = NA))
  },
  # Exon20_Ins
  Exon20_Ins = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["Exon20_Ins"], col = NA))
  },
  # T790M
  T790M = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["T790M"], col = NA))
  },
  # G12C
  G12C = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["G12C"], col = NA))
  },
  # V600E
  V600E = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["V600E"], col = NA))
  },
  # Exon20
  Exon20 = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["Exon20"], col = NA))
  },
  # no_hot
  no_hot = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = vc_cols["no_hot"], col = NA))
  },
  # Linear
  Linear = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["Linear"], col = NA))
  },
  # Complex
  Complex = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["Complex"], col = NA))
  },
  # BFB
  BFB = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["BFB"], col = NA))
  },
  # ecDNA
  ecDNA = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["ecDNA"], col = NA))
  },
  # sv_del
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = vc_cols["DEL"], col = NA))
  },
  # sv_dup
  DUP = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.6, 
              gp = gpar(fill = vc_cols["DUP"], col = NA))
  },
  # sv_inv
  INV = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = vc_cols["INV"], col = NA))
  },
  # CTX
  CTX = function(x, y, w, h) {
    grid.rect(x, y, w*0.4, h*0.4, 
              gp = gpar(fill = vc_cols["CTX"], col = NA))
  },
  # INS
  INS = function(x, y, w, h) {
    grid.rect(x, y, w*0.4, h*0.4, 
              gp = gpar(fill = vc_cols["INS"], col = NA))
  }
)

# correlation heatmap reorder
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

### phenotype colors ------
color_phe <- c('deepskyblue2', '#FFA500', "deepskyblue2", "deepskyblue2", 'plum3', "firebrick3", 'green4', 'gold3', 'plum3', 'gold3', 'plum3',
               'green4', 'gold3', 'plum3', 'darkred', 'yellow2', 'peru', 'dodgerblue4', 'aliceblue', 'mediumorchid3', 'darkcyan', 'lightskyblue1', #1 (data set)
               'darkseagreen3', 'darkseagreen3', "deepskyblue2", 'darkseagreen3', "deepskyblue2", "firebrick3", "firebrick3", "firebrick3", "plum3", "firebrick3", #2 tissue type
               "#F7FBFF", "#DEEBF7", "#C6DBEF", "#C6DBEF", "#C6DBEF", "#9ECAE1", "#9ECAE1", '#FFA500', '#FFA500', '#FFA500', # 3 (8 invasive type)
               brewer.pal(9, "Blues")[2 : 4], '#FFA500', '#FFA500',  #4 (5 invasive type of tumor)
               brewer.pal(9, "Greens")[2 : 4], 'darkseagreen3', 'darkseagreen3', #5 (5 invasive type of normal)
               brewer.pal(9, "Blues")[5], brewer.pal(9, "Blues")[5 : 9], #6 (6 p stage I type)
               brewer.pal(9, "Blues")[5], brewer.pal(9, "Blues")[5 : 9], #7 (6 p stage I type of tumor)
               brewer.pal(9, "Greens")[5], brewer.pal(9, "Greens")[5 : 9], #8 (6 p stage I type of normal)
               brewer.pal(7, "Oranges")[2 : 7], #9 (p stage II_III type)
               brewer.pal(7, "Oranges")[2 : 7], #10 (p stage II_III type of tumor)
               brewer.pal(7, "RdYlGn")[2 : 7], #11 (p stage II_III type of normal)
               brewer.pal(8, "Purples")[4 : 8], '#6BAED6', '#FD8D3C', #12 (p stage)
               "#6BAED6", "#6BAED6", "#6BAED6", "#6BAED6", "#9E9AC8", "#9E9AC8", "#6BAED6", "#9E9AC8", "#FD8D3C", #13 (GGO)
               brewer.pal(9, "Blues")[3 : 9], brewer.pal(9, "Oranges")[3 : 7],  #14 (GGO + invasive stage)
               brewer.pal(9, "Blues")[3], brewer.pal(9, "Purples")[3], brewer.pal(9, "Oranges")[3], brewer.pal(9, "Blues")[5], brewer.pal(9, "Purples")[5], brewer.pal(9, "Oranges")[5],  #15 (AIS/MIA/INV + GGO)
               '#C6DBEF', '#DADAEB', '#FDD0A2', '#6BAED6', '#9E9AC8', '#FD8D3C', '#6BAED6', '#9E9AC8', '#FD8D3C',  #16 (AIS/MIA/INV + GGO)
               brewer.pal(9, "Blues")[3 : 9], brewer.pal(9, "Oranges")[3 : 7], #17 (stage + GGO)
               rep('darkseagreen3', 4), '#C6DBEF', "#6BAED6", brewer.pal(9, "Purples")[5], brewer.pal(9, "Oranges")[5], #18 (tissue type, p stage and GGO)
               'pink', 'deepskyblue2', 'pink', 'deepskyblue2', #19 (gender)
               'lightskyblue1', '#FFA500', '#FFA500', #20 (smoking)
               'salmon', 'darkred', 'lightgrey', #21 (os status)
               'lightgrey', 'black', #22 (mutation status)
               '#cc9e16', 'firebrick3', '#3737f0', '#00CED1', #23 MFP classification
               'firebrick3',  'firebrick3', 'deepskyblue2', 'deepskyblue2', '#FD8D3C', '#FDD0A2',  #24 TP53 mut or wt 
               '#6BAED6', '#FD8D3C', '#6BAED6', '#FD8D3C', '#9E9AC8',  '#FDD0A2', 'gold3', '#2171B5', 'grey', '#6BAED6', '#FDD0A2', '#9E9AC8',  'gold3', '#FD8D3C', 'gold3', '#FD8D3C',  # 25 pathology
               brewer.pal(9, "Purples")[5], brewer.pal(9, "Oranges")[5], # 26 invasive non-solid and solid
               'firebrick3', 'dodgerblue4', 'lightgrey', 'firebrick3', 'dodgerblue4', 'lightgrey', # 27 gain, loss, none
               'firebrick3', 'dodgerblue4', # 28 WGD, Yes, No
               '#FD8D3C', 'dodgerblue4', 'darkseagreen3', '#D94801', # 29 Genes, COL11A1, THBS2
               '#FD8D3C', 'dodgerblue4', # 30 clinical_benefit
               'darkseagreen3', 'dodgerblue4', 'gold3', '#FD8D3C', 'firebrick3', # 31 overall_response 'SR/CR', 'PR', 'SD', 'NE', 'PD'
               '#C6DBEF', 'deepskyblue2', '#9ECAE1', '#2171B5', '#08306B', 'deepskyblue2', '#9ECAE1', '#2171B5', '#D94801', 'yellow2', '#D94801', '#D94801', 'firebrick3', 'firebrick3', '#D94801', 'firebrick3', 'grey', # 32 complex class
               '#C6DBEF', 'deepskyblue2', '#2171B5', '#D94801', 'yellow2', '#D94801', 'firebrick3', 'plum3', '#4A1486', 'firebrick3', # 33 complex class
               '#C6DBEF', 'deepskyblue2', '#2171B5', 'dodgerblue4', '#D94801', 'yellow2', '#4A1486', 'firebrick3',  # 33 complex class
               "#4DAF4A", "#984EA3", "#377EB8", "#E41A1C", # 34 clone type
               "firebrick3", "firebrick3", "#FDAE6B", "#08519C", "deepskyblue2", "deepskyblue2", # 35 relapse
               'lightskyblue1', "#9ECAE1", "#377EB8", 'dodgerblue4', # 36 age
               '#4A1486', "deeppink4", "deeppink4", '#FFA500', 'plum3', 'plum3', 'dodgerblue4', # 37 relapse-site
               'lightskyblue1', 'lightskyblue1', 'pink', 'pink', # 37 relapse-site
               'peru', 'darkred', 'green4', 'lightgrey', 'grey40', # 37 relapse-site
               '#2171B5', '#FFA500', "firebrick3", 'firebrick3', 'firebrick3', 'firebrick3', 'firebrick3',  'firebrick3', 'firebrick3', 'firebrick3', 'firebrick3', # 38 relapse-p
               "#377EB8", "#E41A1C", "#377EB8", "#E41A1C", "#377EB8", "#E41A1C", "#377EB8", "#E41A1C", # 39 relapse-histology
               'deepskyblue2', 'gold3', '#FD8D3C', 'firebrick3', '#4A1486', 'plum3', 'yellow2', '#00CED1', 'darkcyan', 'pink', 'black', 'red', 'green3', 'deepskyblue2', 'gold3', '#FD8D3C', 'firebrick3', '#4A1486', # 40 lacms 
               '#FD8D3C', '#4A1486', # 41 pan pos and neg 
               'deepskyblue2', 'gold3', '#FD8D3C', 'firebrick3', '#4A1486', '#00CED1', '#FDD0A2', # 42 deconvolution methods
               'deepskyblue2', "#FFA500", "#FFA500", "#FFA500", "firebrick3", 'green4', 'gold3', 'plum3', 'lightgrey', # 43 therapy
               'deepskyblue2', 'gold3', '#FD8D3C', 'firebrick3', '#4A1486', 'plum3', 'yellow2', '#00CED1', 'darkcyan', 'pink', # 44 SV subtype  
               'deepskyblue2', 'gold3', '#FD8D3C', 'firebrick3', '#4A1486', 'plum3', 'yellow2', '#00CED1', 'darkcyan', 'pink', # 45 CNV subtype  
               'grey', '#33A02C', '#1F78B4', '#FF7F00', '#E31A1C', # 46 amplicon type
               'grey', 'gold3', '#E31A1C', # 47 chromothripsis
               "deeppink4", "orange3", "purple3", "gray33", "forestgreen", # 48 sv type
               '#FFA500', "deepskyblue2", "firebrick3", 'green4', 'gold3', 'plum3', 'darkred', 'yellow2', # 49 sv signature
               '#FD8D3C', 'firebrick3', 'plum3', '#4A1486', # 50 LUAD location
               'deepskyblue2', 'gold3', '#FD8D3C', 'firebrick3', 'firebrick3', 'firebrick3', 'forestgreen', 'deepskyblue2', '#FF7F00', 'plum3', '#9ECAE1', 'darkseagreen3', '#9ECAE1', 'grey', # 51 SNV
               'darkseagreen3', brewer.pal(9, "Blues")[3:9], brewer.pal(9, "Oranges")[2:9], # 52
               '#C6DBEF', 'deepskyblue2', '#9E9AC8', 'gold3', 'gold3', 'deeppink4', '#FD8D3C', '#FD8D3C',  'firebrick3') # 53 press

 
names(color_phe) <- c('300_GIL', '3000_GIL')