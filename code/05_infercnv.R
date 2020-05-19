library(infercnv)
#setwd('/users/whou/aml_feinberg/infercnv')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/infercnv')
infercnv_obj = CreateInfercnvObject(raw_counts_matrix='genematrix.txt',
                                    annotations_file='meta.txt',
                                    delim="\t",
                                    max_cells_per_group = 500,
                                    gene_order_file='gr.txt',
                                    ref_group_names=c("CD34",'GMP')) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.2, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='output_cutoff0.1_sub500', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)



