# -------------------------------------------
# common CpG between training and AML samples
# -------------------------------------------
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/combine/wgbs/filterme.rds') ## filter out CpG sites with NA on samples, num [1:5585327, 1:284]
amlme.cs <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/aml_feinberg/wgbs/processed/hg38/all_cpg_by_14_patients_completecases.rds')
int <- intersect(rownames(d), rownames(amlme.cs)) ## 3474682
amlme.cs <- amlme.cs[int, ]

# ----------------------------------
# highly variable CpG in AML samples
# ----------------------------------
cm <- rowMeans(amlme.cs)
csv <- sqrt(rowMeans(amlme.cs * amlme.cs - cm^2) / (ncol(amlme.cs) - 1) * ncol(amlme.cs))
# > summary(csv)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.05577 0.09528 0.12206 0.16167 0.51887 
amlme.cs <- amlme.cs[csv >= 0.2, ]
# > str(amlme.cs)
#  num [1:647420, 1:14] 0.833 0.857 0.571 0.875 1 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:647420] "chr1_16243" "chr1_135150" "chr1_184163" "chr1_184175" ...
#   ..$ : chr [1:14] "CD34-1" "SU048" "SU070" "SU204" ...
set.seed(12345)
v = sample(rownames(amlme.cs), 1e4)
saveRDS(d[v, ],file='/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/sampleme_match_aml.rds')
