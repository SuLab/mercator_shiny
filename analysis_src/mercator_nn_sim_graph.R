library(RANN)
library(igraph)
library(philentropy)
library(sets)
library(arrangements)
library(parallel)

getSimilarity <- function(x){

    y <- nn.recount[x[1],]
    z <- nn.recount[x[2],]

    union <- length(union(z,y))
    inter <- length(intersect(z,y))

    sim <- inter/union
    if(sim < 1/15){
        sim <- 0
    }

    return(sim)
}

setwd('project_dir')

nn.recount <- read.table('results/recount/clustering/recount_over50_pc3sd_tpm_log_nv100_90th_var_genes_nn_k60.csv',sep=',',stringsAsFactors=F,row.names=1)
nn.recount <- nn.recount+1

print('beginning')

neighbor.jacc <- matrix(1,nrow=nrow(nn.recount),ncol=nrow(nn.recount))

no_cores <- 9
cl <- makeCluster(no_cores,outfile='/gpfs/group/su/lhgioia/map/results/logs/louvain_clustering_tpm_log_k60.log')

clusterExport(cl,'nn.recount')

print(nrow(nn.recount))

ind.combs <- combinations(nrow(nn.recount),2,replace=T)

weights <- parApply(cl,ind.combs,1,getSimilarity)

print('calculated similarities')

stopCluster(cl)

out.dat <- cbind(ind.combs-1,weights)

write.table(out.dat,'results/recount/clustering/nn_sim_weights_recount_over50_pc3sd_tpm_log_nv100_90th_var_genes_k60.tsv',sep='\t',col.names=F,row.names=F)

