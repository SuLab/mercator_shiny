library(jsonlite)
library(parallel)
library(arrangements)
library(SummarizedExperiment)
library(MASS)

manualNBTest <- function(x){

    library(MASS)

    test.dat <- data.frame(GENE=x,cluster=factor(cluster.ass.sub))

    gene.var <- (var(test.dat$GENE) == 0)

    gene.mean <- (mean(test.dat$GENE) <= 1)

    if(gene.var){
        return(c(2,NA,NA))
    }

    if(gene.mean){
        return(c(5,NA,NA))
    }


    first.missing.vals.test <- sum(subset(test.dat,cluster==names(clus.tab)[1])$GENE) == 0 
    second.missing.vals.test <- sum(subset(test.dat,cluster==names(clus.tab)[2])$GENE) == 0 

    if(first.missing.vals.test & second.missing.vals.test){
        return(c(3,NA,NA))
    } else if(first.missing.vals.test){
        test.dat[test.dat$cluster==names(clus.tab)[1],'GENE'] <- 1
    } else if(second.missing.vals.test){
        test.dat[test.dat$cluster==names(clus.tab)[2],'GENE'] <- 1
    }

    nb.coef <- tryCatch({
        nb.test <- suppressWarnings(glm.nb('GENE ~ cluster',test.dat))
        return(c(summary(nb.test)$coef[2,4],summary(nb.test)$coef[2,1],summary(nb.test)$coef[2,3]))
        },
        error = function(e) {
            return(c(4,NA,NA))
        })

    return(nb.coef)

}


setwd('project_dir')

# make raw count matrix
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

print(dim(cnt.mat))

deseq.sizefactors <- readRDS('data/recount/markers/deseq2_size_factors_rescaled_max_poscounts_90th_var_genes.RDS')

cnt.mat <- cnt.mat[,names(deseq.sizefactors)]

cnt.mat <- t(cnt.mat) / deseq.sizefactors

print(dim(cnt.mat))

cluster.ass <- readRDS('results/recount/clustering/leiden_r9e-3_over50_pc3sd_tpm_log_k30_sim_90th_var_genes.RDS')

cnt.mat <- cnt.mat[intersect(rownames(cnt.mat),names(cluster.ass)),]
cluster.ass <- cluster.ass[rownames(cnt.mat)]

no.cores <- 15
cl <- makeCluster(no.cores,outfile='/gpfs/group/su/lhgioia/map/results/logs/leiden_pairwise_markers_r9e-3_poscounts_nv100.log')

cluster.num <- length(unique(cluster.ass))

if(any(is.na(cluster.ass))){
    cluster.num <- cluster.num-1
}

pval.list <- list()

### we split this up into many jobs to speed up computation on our cluster, you can ignore this splitting 

args <- commandArgs(trailingOnly=TRUE)

curr.comb <- as.integer(args[1]) * 500
file.str <- sprintf('%d_%d',curr.comb,curr.comb + 500)

print(curr.comb)
print(cluster.num)

pval.res <- list()

ind.combs <- combinations(cluster.num,2,replace=F)

for(i in (curr.comb+1):(curr.comb+500)){

    comb <- (ind.combs[i,]-1)

    cluster.ass.sub <- cluster.ass[(cluster.ass == comb[1] | cluster.ass == comb[2])]
    cluster.ass.sub <- sort(cluster.ass.sub)

    clus.tab <- table(cluster.ass.sub)

    cluster.dat <- cnt.mat[names(cluster.ass.sub),]

    clusterExport(cl,'cluster.ass.sub')
    clusterExport(cl,'clus.tab')

    results <- parCapply(cl,cluster.dat,manualNBTest)
    results <- matrix(results,nrow=ncol(cluster.dat),byrow=T)
    rownames(results) <- colnames(cnt.mat)
    colnames(results) <- c('pval','fc','test.stat')

    pval.list[[sprintf('%d.%d',comb[1],comb[2])]] <- results

    print(comb)

    if(i %% 100 == 0){
        print(i)
        saveRDS(pval.list,sprintf('results/recount/markers/tpmlogclus_nv100_90th_var_genes_k30_nb_poscounts90th/leiden_pairwise_markers_pc3sd_poscounts_NB_%s_redo.RDS',file.str))
        next
    }
    if(i == nrow(ind.combs)){
        print(i)
        saveRDS(pval.list,sprintf('results/recount/markers/tpmlogclus_nv100_90th_var_genes_k30_nb_poscounts90th/leiden_pairwise_markers_pc3sd_poscounts_NB_%s_redo.RDS',file.str))
        break
    }
}

saveRDS(pval.list,sprintf('results/recount/markers/tpmlogclus_nv100_90th_var_genes_k30_nb_poscounts90th/leiden_pairwise_markers_pc3sd_poscounts_NB_%s_redo.RDS',file.str))

stopCluster(cl)
