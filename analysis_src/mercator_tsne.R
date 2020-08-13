#################################################
### Mercator t-SNE script
### does a sweep of PC's and perplexities to 
### determine optimal paramters
### uses Rtsne package
#################################################


library(Rtsne)

setwd('project_dir')

set.seed(20)

perplexity.vec <- c(30,40,60,80,110,140,175,210,250)

dim.vec <- c(50,100,250,500,750)

recount.data <- readRDS('results/recount/pca/recount_750_dim_noProj_over50_pc3sd_tpm_log_90th_var_genes.RDS')

NA.rows <- read.table('data/recount/metadata/runs_not_in_metasra.csv',sep=',',header=F,stringsAsFactors=F)

not.NAs <- rownames(recount.data)[!rownames(recount.data) %in% NA.rows$V1]

recount.data <- recount.data[not.NAs,]

print(dim(recount.data))

recount.data <- normalize_input(recount.data)

results <- list()

for(perp in perplexity.vec){
    for(dim in dim.vec){
        tsne.out <- Rtsne(recount.data[,1:dim],perplexity=perp,check_duplicates=FALSE,pca=FALSE,Y_init = recount.data[,1:2],eta=2875)
        rownames(tsne.out$Y) <- rownames(recount.data)
        results[[sprintf('%d.%d',perp,dim)]] <- tsne.out$Y

        saveRDS(results,'results/recount/tsne/recount_tsne_list_pca_pcinit_eta2875_over50_pc3sd_tpm_log_90th_var_genes.RDS')
    }
}
