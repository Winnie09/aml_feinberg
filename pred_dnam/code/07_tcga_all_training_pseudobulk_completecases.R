# -------------------------------------------
# common CpG between training and AML samples
# -------------------------------------------
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/tcga_450k/proc/hg38/combine/me.rds') ## [1:394363, 1:8578]
d <- d[!is.na(rowMeans(d)), ] ## [1:281329, 1:8578] 
loc <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/proc/hg38/loc/hg38.rds')
rownames(d) <- loc[rownames(d)]
amlme.cs <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/wgbs/processed/hg38/all_cpg_by_14_patients_completecases.rds') ## [1:36559571, 1:14]
int <- intersect(rownames(d), rownames(amlme.cs))  ## 106981
amlme.cs <- amlme.cs[int, ]
# ----------------------------------
# highly variable CpG in AML samples
# ----------------------------------
cm <- rowMeans(amlme.cs)
csv <- sqrt(rowMeans(amlme.cs * amlme.cs - cm^2) / (ncol(amlme.cs) - 1) * ncol(amlme.cs))
amlme.cs <- amlme.cs[csv >= 0.2, ]
set.seed(12345)
v = sample(rownames(amlme.cs), 1e4)
meth <- d[v, ]
## ---------------------------
## load, process training data
## ---------------------------
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
source('/home-4//whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/tcga_450k/proc/hg38/combine/ge.rds') ## [1:58560, 1:8578]
project <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/tcga_450k/proc/hg38/combine/project.rds')
meth = meth[complete.cases(meth), ]

# # select LAML samples as training only <<
# meth <- meth[, names(project)[project=='TCGA-LAML']]

expr <- expr[, colnames(meth)]
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ds <- sub(':.*','',colnames(meth)) ##

# avoid 0 and 1 to enter logit function
meth[meth==0] <- min(meth[meth>0])
meth[meth==1] <- max(meth[meth<1])
meth <- log(meth/(1-meth))

# Filtering low expressed genes and quantile normalization
expr <- expr[rowMeans(expr > 0) >= 0.01,]

# Filtering low variable genes
m <- rowMeans(expr)
s <- sqrt((rowMeans(expr*expr) - m^2) * ncol(expr)/(ncol(expr)-1))
mod <- loess(s~m)
expr <- expr[resid(mod) > 0,]
  
## ---------------------------
## load, process testing data 
## ---------------------------
seur <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/scrna/AML_integrated_scRNA_seurat.rds')
predexpr <- as.matrix(seur@assays$RNA@counts)
meta <- seur@meta.data$sample
predexpr <- sapply(unique(meta),function(i) rowSums(predexpr[,meta==i]))

gidgn <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg38_geneid_genename_genelength.rds')
rownames(expr) <- gidgn[match(rownames(expr), gsub('\\..*', '', gidgn$geneid)), 2]
expr <- expr[!duplicated(rownames(expr)), ]

## match genes
intgene <- intersect(rownames(predexpr), rownames(expr))
predexpr <- predexpr[intgene, ]
expr <- expr[intgene, ]

## train model
m <- trainmodel(expr, meth, log.transform = F, filter.low.express.gene = F, filter.low_var.gene = F)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pseudobulk/tcga_training/all_trainmodel_completecase.rds')

## predict
pred <- predict(predexpr, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pseudobulk/tcga_training/all_pred_completecase.rds')

## evaluate
amlme <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/wgbs/processed/hg38/all_cpg_by_14_patients_completecases.rds')
amlme <- amlme[rownames(meth), intersect(colnames(pred), colnames(amlme))]
pred <- pred[, colnames(amlme)]
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)
  }
}

cv1 <- corfunc(pred[complete.cases(amlme),], amlme[complete.cases(amlme), ])
print(summary(cv1))
cv2 <- corfunc(t(pred[complete.cases(amlme), ]),t(amlme[complete.cases(amlme), ]))
print(summary(cv2))

saveRDS(list(crosssample = cv1, crosssite = cv2), '/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pseudobulk/tcga_training/all_perf_completecase.rds')
