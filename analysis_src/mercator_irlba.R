library(irlba)
library(bigmemory)
library(bigalgebra)
library(Rcpp)

matmul <- function(A,B,transpose=FALSE){

    if(is.null(dim(B))) B <- cbind(B)

    if(is.null(dim(A))) A <- rbind(A)

    if(transpose)
        return(cbind((t(B) %*% A)[]))

    ## print(typeof(A))
    ## print(typeof(B))

    ## print(dim(A))
    ## print(dim(B))

    cbind((A %*% B)[])

}

tcga.tpm <- as.matrix(readRDS('data/recount/tcga/tpm_tcga.RDS'))
gtex.tpm <- as.matrix(readRDS('data/recount/gtex/tpm_gtex.RDS'))

## tpm.mat <- rbind(tcga.tpm,gtex.tpm)

tpm.mat <- matrix(0,nrow=70603,ncol=58037)

cnt <- 1

tpm.mat[cnt:nrow(tcga.tpm),] <- tcga.tpm
cnt <- cnt + nrow(tcga.tpm)

tpm.rownames <- rownames(tcga.tpm)

tpm.mat[cnt:(cnt + nrow(gtex.tpm)-1),] <- gtex.tpm
cnt <- cnt + nrow(gtex.tpm)

tpm.rownames <- c(tpm.rownames,rownames(gtex.tpm))

## over50.projs <- readRDS('data/recount/recount_entries_over50_noSingle.RDS')

over50.projs <- read.table('data/recount/recount_over50_labelled.csv',sep=',',stringsAsFactors=F)
over50.projs <- subset(over50.projs,V2=='RNA-seq')


## for(file.id in list.files('data/recount/project_cnts')){
for(file.id in over50.projs$V1){

    file.tpm <- as.matrix(readRDS(sprintf('data/recount/project_cnts/%s/gene_counts_tpm.RDS',file.id)))

    tpm.mat[cnt:(cnt + nrow(file.tpm) - 1),] <- file.tpm
    cnt <- cnt + nrow(file.tpm)

    tpm.rownames <- c(tpm.rownames,rownames(file.tpm))
    ## tpm.mat <- rbind(tpm.mat,file.tpm)

}

print(cnt)
print(dim(tpm.mat))

tpm.mat <- tpm.mat[1:length(tpm.rownames),]

print(dim(tpm.mat))

rownames(tpm.mat) <- tpm.rownames
colnames(tpm.mat) <- colnames(tcga.tpm)
## genes.used <- readRDS('data/recount/new_old_conversion_recount.RDS')

## genes.used <- readRDS('data/recount/recount_protein_coding.RDS')

## colnames(tpm.mat) <- gsub('[.].*$','',colnames(tcga.tpm))
## tpm.mat <- tpm.mat[,genes.used$V3]

NA.rows <- read.table('data/recount/metadata/runs_not_in_metasra.csv',sep=',',header=F,stringsAsFactors=F)

not.NAs <- tpm.rownames[!tpm.rownames %in% NA.rows$V1]

tpm.mat <- tpm.mat[not.NAs,]

print(dim(tpm.mat))

pc3sd.dat <- readRDS('data/recount/recount_over50_bulkOnly_dropout_filtered_pc_3sd.RDS')

pc3sd.samps <- intersect(not.NAs, rownames(pc3sd.dat))
tpm.mat <- tpm.mat[pc3sd.samps,]


var.genes <- readRDS('data/recount/gene_lists/recount_over50_pc3sd_tpm_log_90th_perc_var_genes.RDS')
tpm.mat <- tpm.mat[,var.genes]

print(dim(tpm.mat))

facs <- 1e6/rowSums(tpm.mat)

tpm.mat <- tpm.mat * facs

print(head(rowSums(tpm.mat)))

tpm.mat <- log(tpm.mat+1)

## var.genes <- readRDS('data/recount/gene_lists/recount_over50_pc3sd_tpm_log_90th_perc_var_genes.RDS')
## tpm.mat <- tpm.mat[,var.genes]

## print('making big matrix')

## big.data <- as.big.matrix(tpm.mat,backingfile='all_tpm_noProj.bin',descriptorfile='all_tpm_noProj.desc',backingpath='/gpfs/group/su/lhgioia/map/results/recount/pca/tmp')

## for(i in 1:ncol(big.data)){
##     big.data[,i] - big.data[,i] - mean(big.data[,i])
## }

col.means <- c()
col.vars <- c()

for(i in 1:ncol(tpm.mat)){
    col.means <- c(col.means,mean(tpm.mat[,i]))
    col.vars <- c(col.vars,sd(tpm.mat[,i]))
    
    tpm.mat[,i] <- tpm.mat[,i] - mean(tpm.mat[,i])
    tpm.mat[,i] <- tpm.mat[,i] / sd(tpm.mat[,i])
}

param.df <- data.frame(means=col.means,vars=col.vars)
## rownames(param.df) <- colnames(tcga.tpm)
rownames(param.df) <- colnames(tpm.mat)
saveRDS(param.df,'results/recount/pca/tmp/recount_over50_pc3sd_colparams_tpm_log_90th_var_genes.RDS')

## tpm.mat <- tpm.mat[,order(col.vars,decreasing=T)[1:10000]]

## tpm.mat <- tpm.mat[,col.vars > 0]

## names(col.means) <- colnames(tpm.mat)
## saveRDS(col.means,'/gpfs/group/su/lhgioia/map/results/recount/pca/tmp/recount_over50_noSingle_centers.RDS')

## print('attaching big matrix')

## big.data <- attach.big.matrix('/gpfs/group/su/lhgioia/map/results/recount/pca/tmp/all_tpm_noProj.desc')

print ('beginning pca')

gene.pca <- irlba(tpm.mat,nv=750,nu=0,mult=matmul)

saveRDS(gene.pca,'results/recount/pca/recount_all_big_irlba_750_over50_pc3sd_tpm_log_90th_var_genes.RDS')

## gene.pca <- readRDS('/gpfs/group/su/lhgioia/map/results/recount/pca/recount_all_big_irlba.RDS')

rot.gene.dat <- as.matrix(tpm.mat %*% gene.pca$v)
rownames(rot.gene.dat) <- rownames(tpm.mat)

saveRDS(rot.gene.dat,'recount/pca/recount_750_dim_noProj_over50_pc3sd_tpm_log_90th_var_genes.RDS')
