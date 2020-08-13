###############################################################
### Simply do one-vs-all calculations for marker information
###############################################################

library(jsonlite)
library(parallel)
library(arrangements)
library(SummarizedExperiment)

manualTTest <- function(x){

    first.clust.dat <- dist.info[[x]]
    second.clust.dat <- dist.comp.info[[x]]

    diff.stat <- (first.clust.dat$means - second.clust.dat$means) / sqrt(first.clust.dat$vars/first.clust.dat$n + second.clust.dat$vars / second.clust.dat$n)

    diff.dof <- (first.clust.dat$vars / first.clust.dat$n + second.clust.dat$vars / second.clust.dat$n)^2 / (first.clust.dat$vars^2 / (first.clust.dat$n^2 * (first.clust.dat$n-1)) + second.clust.dat$vars^2 / (second.clust.dat$n^2 * (second.clust.dat$n-1)))

    p.vals <- sapply(1:length(diff.stat),function(x) pt(-1*abs(diff.stat[x]),diff.dof[x])*2)

    zero.vars <- (first.clust.dat$vars == 0) & (second.clust.dat$vars == 0)
    p.vals[zero.vars] <- 2

    total.diff <- first.clust.dat$means - second.clust.dat$means

    res.df <- data.frame(pvalue = p.vals,
                         fold.change = total.diff)

    rownames(res.df) <- names(first.clust.dat$means)

    print(x)

    return(res.df)

}


print('beginning')

setwd('/gpfs/group/su/lhgioia/map/')

load('data/recount/gtex/rse_gene_gtex.Rdata')
gtex.dat <- assays(rse_gene)$counts

load('data/recount/tcga/rse_gene_tcga.Rdata')
tcga.dat <- assays(rse_gene)$counts

cnt.mat <- matrix(0,nrow=58037,ncol=70603)

cnt <- 1

cnt.mat[,cnt:ncol(tcga.dat)] <- tcga.dat
cnt <- cnt + ncol(tcga.dat)

cnt.colnames <- colnames(tcga.dat)

cnt.mat[,cnt:(cnt + ncol(gtex.dat)-1)] <- gtex.dat
cnt <- cnt + ncol(gtex.dat)

cnt.colnames <- c(cnt.colnames,colnames(gtex.dat))

over50.projs <- read.table('data/recount/recount_over50_labelled.csv',sep=',',stringsAsFactors=F)
over50.projs <- subset(over50.projs,V2=='RNA-seq')

for(file.id in over50.projs$V1){

    load(sprintf('data/recount/project_cnts/%s/rse_gene.Rdata',file.id))
    proj.dat <- assays(rse_gene)$counts
    
    cnt.mat[,cnt:(cnt + ncol(proj.dat) - 1)] <- proj.dat
    cnt <- cnt + ncol(proj.dat)

    cnt.colnames <- c(cnt.colnames,colnames(proj.dat))

}

print(cnt)
print(dim(cnt.mat))

cnt.mat <- cnt.mat[,1:length(cnt.colnames)]

print(dim(cnt.mat))

colnames(cnt.mat) <- cnt.colnames

NA.rows <- read.table('data/recount/metadata/runs_not_in_metasra.csv',sep=',',header=F,stringsAsFactors=F)

not.NAs <- cnt.colnames[!cnt.colnames %in% NA.rows$V1]
cnt.mat <- cnt.mat[,not.NAs]

print(dim(cnt.mat))

pc3sd.dat <- readRDS('data/recount/recount_over50_bulkOnly_dropout_filtered_pc_3sd.RDS')
pc3sd.samps <- intersect(not.NAs, rownames(pc3sd.dat))
cnt.mat <- cnt.mat[,pc3sd.samps]

print(dim(cnt.mat))

rownames(cnt.mat) <- rownames(proj.dat)

for(i in 1:ncol(cnt.mat)){
    if(max(cnt.mat[,i]) > .Machine$integer.max){
        cnt.mat[,i] <- cnt.mat[,i] * (.Machine$integer.max) / (max(cnt.mat[,i]))
    }
}

gene.vars <- apply(cnt.mat,1,var)
cnt.mat <- cnt.mat[gene.vars > 0,]

deseq.sizefactors <- readRDS('data/recount/markers/deseq2_size_factors_rescaled_max_poscounts_90th_var_genes.RDS')

cnt.mat <- cnt.mat[,names(deseq.sizefactors)]

cnt.mat <- t(cnt.mat) / deseq.sizefactors

print(dim(cnt.mat))
 

cluster.ass <- readRDS('results/recount/clustering/leiden_r9e-3_over50_pc3sd_tpm_log_k30_sim_90th_var_genes.RDS')
cluster.ass <- cluster.ass[rownames(cnt.mat)]

cluster.num <- length(unique(cluster.ass))

if(any(is.na(cluster.ass))){
    cluster.num <- cluster.num-1
}

print(cluster.num)

dist.comp.info <- list()

dist.info <- list()

for(i in 1:cluster.num){

    print(sprintf('calculating info for cluster %d',i))

    clust.samps <- intersect(rownames(cnt.mat),names(which(cluster.ass==i-1)))
    clust.dat <- cnt.mat[clust.samps,,drop=F]

    clust.means <- apply(clust.dat,2,mean)
    clust.vars <- apply(clust.dat,2,var)
    dropout.dat <- apply(clust.dat,2,function(x) sum(x > 0))

    dist.info[[i]] <- list(means=clust.means,vars=clust.vars,dropouts=dropout.dat,n=nrow(clust.dat))

    comp.clust.samps <- intersect(rownames(cnt.mat),names(which(cluster.ass!=i-1)))
    comp.clust.dat <- cnt.mat[comp.clust.samps,,drop=F]

    comp.clust.means <- apply(comp.clust.dat,2,mean)

    clust.means <- log2(clust.means+1)
    clust.means[is.infinite(clust.means)] <- 0
    
    comp.clust.means <- log2(comp.clust.means+1)
    comp.clust.means[is.infinite(comp.clust.means)] <- 0

    entry <- clust.means - comp.clust.means

    dist.comp.info[[i]] <- entry

}

saveRDS(dist.comp.info,'data/recount/markers/leiden_tpmlog_clus_pca_over50_pc3sd_poscounts_comp_dist_info.RDS')
    
