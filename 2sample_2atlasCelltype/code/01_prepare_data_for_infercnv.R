library(data.table)
gtf <- fread('/dcl02/hongkai/data/whou/resource/gtf/grch38.gtf',data.table=F)  ## grch38 = hg38（13161, more genes intersect with hg19）
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\";.*','',sub('.*; gene_name \"','',gtf[,9]))

## prepare gene by cell expression matrix
d1 <- readRDS('/dcl02/hongkai/data/wzhou14/Andy_lab/AML_project/scRNA_process/CD34-1_scRNA_seurat.rds')
cnt1 <- as.matrix(d1@assays$RNA@counts)
cnt1 = cnt1[!duplicated(rownames(cnt1)), ]
colnames(cnt1) <- paste0('CD34:',colnames(cnt1))
d2 <- readRDS('/dcl02/hongkai/data/wzhou14/Andy_lab/AML_project/scRNA_process/GMP-1_scRNA_seurat.rds')
cnt2 <- as.matrix(d2@assays$RNA@counts)
cnt2 = cnt2[!duplicated(rownames(cnt2)), ]
colnames(cnt2) <- paste0('GMP:',colnames(cnt2))
intgene <- intersect(rownames(cnt1), rownames(cnt2))
cnt <- cbind(cnt1[intgene,], cnt2[intgene,])

d3 <- readRDS('/dcl02/hongkai/data/wzhou14/Andy_lab/AML_project/scRNA_process/SU344_scRNA_seurat.rds')
cnt3 <- as.matrix(d3@assays$RNA@counts)
cnt3 = cnt3[!duplicated(rownames(cnt3)), ]
colnames(cnt3) <- paste0('SU344:',colnames(cnt3))
d4 <- readRDS('/dcl02/hongkai/data/wzhou14/Andy_lab/AML_project/scRNA_process/SU462_scRNA_seurat.rds')
cnt4 <- as.matrix(d4@assays$RNA@counts)
cnt4 = cnt4[!duplicated(rownames(cnt4)), ]
colnames(cnt4) <- paste0('SU462:',colnames(cnt4))
intgene <- intersect(rownames(cnt3), rownames(cnt4))
cnt_aml <- cbind(cnt3[intgene,], cnt4[intgene,])

intgene <- intersect(rownames(cnt), rownames(cnt_aml))
cnt <- cbind(cnt[intgene,], cnt_aml[intgene,])

## cell annotation file
meta <- data.frame(cell = colnames(cnt), sample = sub(':.*','', colnames(cnt)), stringsAsFactors = F)
write.table(meta,file='/users/whou/aml_feinberg/infercnv/meta.txt',quote=F,sep='\t',col.names=F,row.names=F)

## gene annotation file
gr <- data.frame(gene = row.names(cnt),gtf[match(row.names(cnt),gn),c(1,4,5)],stringsAsFactors = F)
colnames(gr) = c('gene','chr','start','end')
gr <- gr[!is.na(gr[,4]),]
write.table(gr,file='/users/whou/aml_feinberg/infercnv/gr.txt',quote=F,sep='\t',col.names=F,row.names=F)

## write gene by cell expression matrix
cnt <- cnt[gr[,1], ]
write.table(cnt,file='/users/whou/aml_feinberg/infercnv/genematrix.txt',quote=F,sep='\t')

