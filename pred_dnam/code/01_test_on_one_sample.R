## ---------------------------
## load, process training data
## ---------------------------
library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
source('/home-4//whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
meth <- readRDS('/home-4//whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/sampleme.rds')
expr <- readRDS('/home-4//whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/filterge.rds')
ds <- sub(':.*','',colnames(meth))

# avoid 0 and 1 to enter logit function
meth[meth==0] <- min(meth[meth>0])
meth[meth==1] <- max(meth[meth<1])
meth <- log(meth/(1-meth))

# Filtering low expressed genes and quantile normalization
expr <- expr[rowMeans(expr > 0) >= 0.01,]
qnem <- rowMeans(apply(expr,2,function(i) sort(i)))
gn <- row.names(expr)
expr <- apply(expr,2,function(i) qnem[frank(i,ties.method='min')])
row.names(expr) <- gn

# Filtering low variable genes
m <- rowMeans(expr)
s <- sqrt((rowMeans(expr*expr) - m^2) * ncol(expr)/(ncol(expr)-1))
mod <- loess(s~m)
expr <- expr[resid(mod) > 0,]

## ---------------------------
## load, process testing data 
## ---------------------------
af = list.files('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/scrna')
expraf = af[grepl('seurat', af) & !grepl('AML', af)]
expras <- gsub('_.*', '', expraf)
methaf <- list.files('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/wgbs')
methas <- gsub('_.*', '', methaf)
as <- intersect(expras, methas)
s = as[1]
exprf <- expraf[expras == s]
methf <- methaf[methas == s]

seur <- readRDS(paste0('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/scrna/', exprf))  
predexpr <- as.matrix(seur@assays$RNA@counts)

# predmeth <- read.table(paste0('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/wgbs/', methf))

# normalize by libsize, log-transform
lb <- colSums(predexpr)
lb <- lb/median(lb)
predexpr <- t(t(predexpr)/lb)
predexpr <- log2(predexpr + 1)
# # > summary(rowMeans(predexpr > 0))
# #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # 0.0002395 0.0059880 0.0728144 0.1751167 0.2661078 1.0000000 
# predexpr <- predexpr[rowMeans(predexpr > 0) >= 0.01, ]

gidgn <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg38_geneid_genename_genelength.rds')
rownames(expr) <- gidgn[match(rownames(expr), gsub('\\..*', '', gidgn$geneid)), 2]
expr <- expr[!duplicated(rownames(expr)), ]

## match genes
intgene <- intersect(rownames(predexpr), rownames(expr))
predexpr <- predexpr[intgene, ]
expr <- expr[intgene, ]

## train model
m <- trainmodel(expr, meth, log.transform = F, filter.low.express.gene = F, filter.low_var.gene = F)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/trainmodel.rds')

## predict

## evaluate
############ continue on prediting 
pred <- predict(predexpr, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/aml_feinberg/pred_dnam/result/pred.rds')


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
cv <- corfunc(pred,meth[,ds %in% testds])
print(summary(cv))
cv <- corfunc(t(pred),t(meth[,ds %in% testds]))
print(summary(cv))
# }

