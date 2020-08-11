tcga_perf <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pseudobulk/tcga_training/all_perf_completecase.rds')

tcga_laml_perf <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pseudobulk/tcga_training/laml_perf.rds')

wgbs_perf <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pseudobulk/wgbs_training/perf.rds')

tcga_laml_perf_traingcv0.2 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pseudobulk/tcga_training/all_perf_completecase_trainingcv0.2.rds')


len = length(tcga_laml_perf[[1]])
library(reshape2)
df1 <- data.frame(tcga = tcga_perf[[1]][1:len], tcga_laml = tcga_laml_perf[[1]], wgbs = wgbs_perf[[1]][1:len], tcga_laml_traingcv0.2 = tcga_laml_perf_traingcv0.2[[1]][1:len])
pd1 <- melt(df1)
df2 <- data.frame(tcga = tcga_perf[[2]], tcga_laml = tcga_laml_perf[[2]], wgbs = wgbs_perf[[2]], tcga_laml_traingcv0.2 = tcga_laml_perf_traingcv0.2[[2]])
pd2 <- melt(df2)


library(ggplot2)
library(RColorBrewer)
p1 <- ggplot(pd1, aes(x = variable, y = value, fill = variable)) +
  geom_jitter(shape='.', position=position_jitter(0.2), alpha = 0.5, size = 0.5) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('') + 
  ylab('Cross-sample Spearman Correlation') +
  theme(axis.text.x = element_text(size = 10, color='black'))
  



p2 <- ggplot(pd2, aes(x = variable, y = value, fill = variable)) +
  geom_jitter(shape='.', position=position_jitter(0.2), alpha = 1, size = 10) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab('') + 
  ylab('Cross-site Spearman Correlation') +
  theme(axis.text.x = element_text(size = 10, color='black'))
  
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/plot/plot/perf.pdf',width=5,height=6)
gridExtra::grid.arrange(p1,p2,nrow=2)
dev.off()
