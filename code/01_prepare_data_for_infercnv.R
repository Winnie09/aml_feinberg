library(data.table)
gtf <- fread('/dcl02/hongkai/data/whou/resource/gtf/grch38.gtf',data.table=F)  ## grch38 = hg38（13161, more genes intersect with hg19）
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\";.*','',sub('.*; gene_name \"','',gtf[,9]))

## prepare reference data: gene by cell count matrix
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

###########
# ddir = '/users/whou/aml_feinberg/data/'
ddir <- '/dcl02/hongkai/data/wzhou14/Andy_lab/AML_project/scRNA_process/'
af <- list.files(ddir)
af = af[grepl('_seurat.rds', af)]
af = af[af != 'AML_integrated_scRNA_seurat.rds']
cntlist <- list()

for (f in af){
  s = sub('_.*','',f)
  print(s)
  d <- readRDS(paste0(ddir, f))@assays$RNA@counts
  d <- as.matrix(d)
  cnt = d[!duplicated(rownames(d)), ]
  colnames(cnt) <- paste0(s,':', colnames(cnt))
  print(str(cnt))
  cntlist[[s]] <- cnt
}


intgene <- names(which(table(unlist(lapply(cntlist, rownames)))==length(cntlist)))

cntlist <- lapply(cntlist, function(i) i[intgene, ])
cnt_aml <- do.call(cbind, cntlist)

############# not yet
intgene <- intersect(rownames(cnt), rownames(cnt_aml))
cnt_all <- cbind(cnt[intgene,], cnt_aml[intgene,])

## cell annotation file
meta <- data.frame(cell = colnames(cnt_all), sample = sub(':.*','', colnames(cnt_all)), stringsAsFactors = F)
write.table(meta,file='/users/whou/aml_feinberg/infercnv/meta.txt',quote=F,sep='\t',col.names=F,row.names=F)

## gene annotation file
gr <- data.frame(gene = row.names(cnt_all),gtf[match(row.names(cnt_all),gn),c(1,4,5)],stringsAsFactors = F)
colnames(gr) = c('gene','chr','start','end')
gr <- gr[!is.na(gr[,4]),]
write.table(gr,file='/users/whou/aml_feinberg/infercnv/gr.txt',quote=F,sep='\t',col.names=F,row.names=F)

## write gene by cell expression matrix
cnt_all <- cnt_all[gr[,1], ]
#write.table(cnt_all,file='/users/whou/aml_feinberg/infercnv/genematrix.txt',quote=F,sep='\t')
df = as.data.frame(cnt_all, stringsAsFactors=F)
fwrite(df, file = '/users/whou/aml_feinberg/infercnv/genematrix.txt', quote = F, sep = '\t', row.names = T)
