library(infercnv)
setwd('/users/whou/aml_feinberg/infercnv')
infercnv_obj = CreateInfercnvObject(raw_counts_matrix='genematrix.txt',
                                    annotations_file='meta.txt',
                                    delim="\t",
                                    gene_order_file='gr.txt',
                                    ref_group_names=c("CD34",'GMP')) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='output', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)



