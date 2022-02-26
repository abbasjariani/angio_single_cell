setwd("~/myProjects/211211_system_biology_baghiatallah/angio_genesis/")
library(reshape2)

library(pheatmap)
library(stringr)
library(Seurat)
library(ggplot2)
library(data.table)
#reticulate::py_install(packages = 'umap-learn')
#download clustering files from
#https://singlecell.broadinstitute.org/single_cell/study/SCP1415/cryopreservation-of-human-cancers-conserves-tumour-heterogeneity-for-single-cell-multi-omics-analysis?genes=PTEN#study-download

#check interaction on hippie interaction


###########################################
clust_data_1 <-
  read.csv("CID4471_clusteringfile.txt", sep = '\t')
#breast 1
clust_data_2 <-
  read.csv("CID44971_clusteringfile.txt", sep = '\t')
clust_data_3 <-
  read.csv("CID4513_clusteringfile.txt", sep = '\t')
clust_data_4 <-
  read.csv("PID17267_clusteringfile.txt", sep = '\t')
clust_data_5 <-
  read.csv("SCC180161_clusteringfile.txt", sep = '\t')

clust_data <-
  rbind(clust_data_1, clust_data_2, clust_data_3, clust_data_4, clust_data_5)


#clust_data$study <-
#  str_split_fixed(clust_data$NAME, '_', n= 2)[,1]

#ggplot(clust_data, aes(X, Y, color = interaction(study, ClusterID )))+
#  geom_point()

##################################3
df_angio_genes <-
  read.csv('angiogenesis2.csv')
gene_vec <- df_angio_genes$Gene
gene_vec <- gsub(' ', '',gene_vec )
gene_vec <-
  gene_vec[gene_vec!='']
###########################################################3
################################################################

data_Wu <-
  Read10X(data.dir = "/home/abbas/myProjects/211211_system_biology_baghiatallah/angio_genesis/wu_sc_rnaseq_data/",
          gene.column = 1)


length(intersect(gene_vec, rownames(data_Wu)))
setdiff(gene_vec, rownames(data_Wu))

write.csv(rownames(data_Wu), 'all_genes.csv')
#genes in the angio genesis list that are not present in 10x
#"PANTR1"  "OR10J5"  "MIR503"  "MIRN378" "MIR126"  "MINAR1"  "MIR92A1" "NAXE"    "CNMD"    "MIR20B"  "CYP4Z2P" "WARS1"  
#[13] "MIR132"  "MIR30E"  "MIR17HG" "PI5"     "EVR1"    "CCN1"    "FOXO1A"  "FOXO3A"  ""        "MACIR"   "MIR558"  "MIR30B" 
#[25] "MIR9-3"  "RNASE9"  "MIR19B1" "GDF2"    "SARS1"   "MIR34A"  "YARS1"   "MIR130A"
#some are non-coding RNA

rownames(data_Wu)[grepl('^OR10J', rownames(data_Wu))]
##################################################

data_Wu_obj <- CreateSeuratObject(counts = data_Wu, project = "wu", min.cells = 3, min.features = 200)

data_Wu_obj@meta.data$NAME <- rownames(data_Wu_obj@meta.data)

data_Wu_obj@meta.data <- merge(data_Wu_obj@meta.data, clust_data[,c('NAME', 'CellType', 'X', 'Y')] )


Idents(object = data_Wu_obj) <- "CellType"

head(
  data_Wu_obj@meta.data )
colnames(data_Wu_obj)

#diff_test_res <-
#  FindMarkers(data_Wu_obj)

for(cur_cell_type in c('B-cells',
                       'cDCs', 
                       'pDCs', 
                       'Myoepithelial', 
                       'T-cells', 
                       'CAFs', 
                       'MAST cells', 
                       'NK cells', 
                       'Endothelial',

                       'Cancer/Epithelial',
                       'Monocytes/Macrophages',
                       'Plasmablasts',
                       'PVL cells')){
  print(cur_cell_type)
  diff_test_res <-
    FindMarkers(data_Wu_obj, ident.1 = cur_cell_type)
  if (cur_cell_type == 'Cancer/Epithelial'){
    write.csv(diff_test_res, paste0('Cancer_Epithelial', '_diff_exp_res.csv'), row.names = T)
  }
  if (cur_cell_type == 'Monocytes/Macrophages'){
    write.csv(diff_test_res, paste0('Monocytes_Macrophages', '_diff_exp_res.csv'), row.names = T)
  }
  else{
    write.csv(diff_test_res, paste0(cur_cell_type, '_diff_exp_res.csv'), row.names = T)
  }
}
######################################



##################################
angio_genes_df <- expand.grid(gene_vec[gene_vec!= ''], c('B-cells',
                                       'cDCs', 
                                       'pDCs', 
                                       'Myoepithelial', 
                                       'T-cells', 
                                       'CAFs', 
                                       'MAST cells', 
                                       'NK cells', 
                                       'Cancer',
                                       'Endothelial',
                                       'Endothelial 1', 
                                       'Endothelial 2',
                                       'Cancer_Epithelial',
                                       'Monocytes_Macrophages',
                                       'Plasmablasts',
                                       'PVL cells'))


angio_genes_df <- as.data.table(angio_genes_df)
colnames(angio_genes_df) <- c('gene', 'cell_type')
angio_genes_df$OE <- 0

for(cur_cell_type in c('B-cells',
  'cDCs', 
  'pDCs', 
  'Myoepithelial', 
  'T-cells', 
  'CAFs', 
  'MAST cells', 
  'NK cells', 
  'Cancer',
  'Endothelial',
  'Endothelial 1', 
  'Endothelial 2',
  'Cancer_Epithelial',
  'Monocytes_Macrophages',
  'Plasmablasts',
  'PVL cells')){
  
  
  
  cur_diff_exp_data <-
    read.csv(paste0("diff_exp_res/", cur_cell_type,'_diff_exp_res.csv'))
  
  cur_diff_exp_data_filt <-
    cur_diff_exp_data[cur_diff_exp_data$p_val_adj< 0.05 & 
                        cur_diff_exp_data$avg_log2FC > 0.2, ]
  
  
  #cur_cell_type <-
  #  gsub(' ', '_' , cur_cell_type)
  
  #angio_genes_df[angio_genes_df$gene%in%cur_diff_exp_data_filt$X, cur_cell_type := 'OE']
  angio_genes_df[angio_genes_df$gene%in%cur_diff_exp_data_filt$X &
                   angio_genes_df$cell_type== cur_cell_type, OE := 1]
  write.csv(intersect(cur_diff_exp_data_filt$X, gene_vec),
            paste0('/home/abbas/myProjects/211211_system_biology_baghiatallah/angio_genesis/angio_genes_OE_cellTypes/',
                   cur_cell_type , '.csv'))
    
}



ggplot(angio_genes_df, aes(x = cell_type, y=gene, fill = OE))+
  geom_tile()
angio_genes_df_cast <-
 dcast(angio_genes_df, gene ~ cell_type , value.var = 'OE')


angio_genes_df_cast_df <- as.data.frame(angio_genes_df_cast)
rownames(angio_genes_df_cast_df) <- angio_genes_df_cast_df$gene
angio_genes_df_cast_df <- angio_genes_df_cast_df[,colnames(angio_genes_df_cast_df)!='gene']

angio_genes_df_cast_df_red <-
  angio_genes_df_cast_df[rowSums(angio_genes_df_cast_df)!=0 , colSums(angio_genes_df_cast_df) >=  0]
angio_genes_df_cast_df_red <-
  angio_genes_df_cast_df_red[,!colnames(angio_genes_df_cast_df_red)%in%c('Endothelial 1','Endothelial 2','Cancer') ]
pheatmap(
  angio_genes_df_cast_df_red, fontsize_row = 2)

#library(ComplexHeatmap)

#Heatmap(angio_genes_df_cast)
######################
genes_group1 <-
  rownames(angio_genes_df_cast_df_red[angio_genes_df_cast_df_red$Endothelial == 1 |
                               angio_genes_df_cast_df_red$`PVL cells` == 1 |
                               angio_genes_df_cast_df_red$Myoepithelial == 1 , ] )
write.csv(genes_group1, 'genes_group1.csv')

genes_group2 <-
  rownames(angio_genes_df_cast_df_red[angio_genes_df_cast_df_red$Monocytes_Macrophages == 1 , ] )
write.csv(genes_group2, 'genes_group2.csv')


genes_group3 <-
  rownames(angio_genes_df_cast_df_red[angio_genes_df_cast_df_red$`NK cells` == 1 |
                                        angio_genes_df_cast_df_red$`MAST cells` == 1 |
                                        angio_genes_df_cast_df_red$Cancer_Epithelial == 1 , ] )
write.csv(genes_group3, 'genes_group3.csv')