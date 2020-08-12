setwd('project_dir')

marker.dat <- list()

for(i in 0:23){

    marker.slice <- readRDS(sprintf('results/recount/markers/tpmlogclus_nb_poscounts/leiden_pairwise_markers_pc3sd_poscounts_NB_%d2%d_redo.RDS',500*i,500*(i+1)))
    
    marker.dat <- c(marker.dat,marker.slice)
    
}


cluster.ass <- readRDS('results/recount/clustering/leiden_r25e-3_pc3sd_tpm_log_k40_sim_nosingles.RDS')


cluster.num <- length(unique(cluster.ass))

if(any(is.na(cluster.ass))){
    cluster.num <- cluster.num-1
}

gene.file <- read.table('data/recount/gene_names_all.txt',stringsAsFactors=F)
gene.ids <- gene.file$V1

cluster.unique.pvals <- list()
cluster.pvals <- list()

cluster.up.fcs <- list()
cluster.min.fcs <- list()

cluster.down.fcs <- list()
cluster.max.fcs <- list()

cluster.min.stats <- list()
cluster.max.stats <- list()

for(i in 1:cluster.num){

    print(i)

    comp.mat <- matrix(1,nrow=length(gene.ids),ncol=cluster.num-1)
    rownames(comp.mat) <- gene.ids

    fc.mat <- matrix(0,nrow=length(gene.ids),ncol=cluster.num-1)
    rownames(fc.mat) <- gene.ids

    stat.mat <- matrix(0,nrow=length(gene.ids),ncol=cluster.num-1)
    rownames(stat.mat) <- gene.ids

    cluster.comps <- names(marker.dat)[grepl(sprintf('^%d[.]',(i-1)),names(marker.dat)) | grepl(sprintf('[.]%d$',(i-1)),names(marker.dat))]

    if(length(cluster.comps) != cluster.num-1){
        print(sprintf('Wrong number of comps for cluster %d',i))
        break
    }

    for(j in 1:length(cluster.comps)){

        comp.vec <- marker.dat[[cluster.comps[j]]]


        comp.mat[rownames(comp.vec),j] <- comp.vec[,'pval']
        if(i <= j){
            fc.mat[rownames(comp.vec),j] <- -1 * comp.vec[,'fc']
            stat.mat[rownames(comp.vec),j] <- -1 * comp.vec[,'test.stat']
        }
       else{
            fc.mat[rownames(comp.vec),j] <- comp.vec[,'fc']
            stat.mat[rownames(comp.vec),j] <- comp.vec[,'test.stat']
            
        }

    }
    
    gene.unique.pvals <- rep(0,nrow(comp.mat))
    gene.pvals <- rep(0,nrow(comp.mat))

    gene.up.fcs <- rep(0,nrow(comp.mat))
    gene.min.fcs <- rep(0,nrow(comp.mat))
    
    gene.down.fcs <- rep(0,nrow(comp.mat))
    gene.max.fcs <- rep(0,nrow(comp.mat))
    
    gene.min.stats <- rep(0,nrow(comp.mat))
    gene.max.stats <- rep(0,nrow(comp.mat))

    for(j in 1:nrow(comp.mat)){

        comp.row <- comp.mat[j,]
        row.sorted <- sort(comp.row,decreasing=TRUE)
        gene.pvals[j] <- row.sorted[8]

        gene.unique.pvals[j] <- row.sorted[1]

        fc.row <- fc.mat[j,]
        fc.sorted <- sort(fc.row,decreasing=FALSE)
        gene.up.fcs[j] <- fc.sorted[8]
        gene.min.fcs[j] <- fc.sorted[1]

        fc.sorted <- sort(fc.row,decreasing=TRUE)
        gene.down.fcs[j] <- fc.sorted[8]
        gene.max.fcs[j] <- fc.sorted[1]

        stat.row <- stat.mat[j,]
        stat.sorted <- sort(stat.row,decreasing=FALSE)
        gene.min.stats[j] <- stat.sorted[1]
        gene.max.stats[j] <- stat.sorted[length(stat.sorted)]

    }

    names(gene.unique.pvals) <- rownames(comp.mat)
    names(gene.pvals) <- rownames(comp.mat)

    names(gene.up.fcs) <- rownames(fc.mat)
    names(gene.min.fcs) <- rownames(fc.mat)
    
    names(gene.down.fcs) <- rownames(fc.mat)
    names(gene.max.fcs) <- rownames(comp.mat)

    names(gene.min.stats) <- rownames(comp.mat)
    names(gene.max.stats) <- rownames(comp.mat)
    
    cluster.pvals[[i]] <- gene.pvals
    cluster.unique.pvals[[i]] <- gene.unique.pvals

    cluster.up.fcs[[i]] <- gene.up.fcs
    cluster.min.fcs[[i]] <- gene.min.fcs

    cluster.down.fcs[[i]] <- gene.down.fcs
    cluster.max.fcs[[i]] <- gene.max.fcs

    cluster.max.stats[[i]] <- gene.max.stats
    cluster.min.stats[[i]] <- gene.min.stats
}



marker.list <- list()

onevsall.info <- readRDS('data/recount/markers/leiden_tpmlog_clus_pca_over50_pc3sd_poscounts_comp_dist_info.RDS')

for(i in 1:length(cluster.pvals)){
    print(i)
    pvals <- p.adjust(cluster.pvals[[i]],method='BH')
    unique.pvals <- p.adjust(cluster.unique.pvals[[i]],method='BH')

    up.fcs <- cluster.up.fcs[[i]]
    min.fcs <- cluster.min.fcs[[i]]
    
    down.fcs <- cluster.down.fcs[[i]]
    max.fcs <- cluster.max.fcs[[i]]
    
    min.stats <- cluster.min.stats[[i]]
    max.stats <- cluster.max.stats[[i]]
    
    markers <- pvals[!is.na(pvals)]
    unique.pvals <- unique.pvals[!is.na(pvals)]

    up.fcs <- up.fcs[!is.na(pvals)]
    min.fcs <- min.fcs[!is.na(pvals)]
    
    down.fcs <- down.fcs[!is.na(pvals)]
    max.fcs <- max.fcs[!is.na(pvals)]
    
    min.stats <- min.stats[!is.na(pvals)]
    max.stats <- max.stats[!is.na(pvals)]

    up.marker.filter <- (markers < 0.01) & (up.fcs > 0)
    down.marker.filter <- (markers < 0.01) & (down.fcs < 0)

    markers <- c(markers[up.marker.filter],markers[down.marker.filter])
    unique.pvals <- c(unique.pvals[up.marker.filter],unique.pvals[down.marker.filter])

    fcs <- c(up.fcs[up.marker.filter],down.fcs[down.marker.filter])

    stats <- c(min.stats[up.marker.filter],max.stats[down.marker.filter])

    unique.pvals <- unique.pvals[unique(names(markers))]
    fcs <- fcs[unique(names(markers))]
    stats <- stats[unique(names(markers))]
    markers <- markers[unique(names(markers))]

    onevsall.entry <- onevsall.info[[i]]
    total.fcs <- onevsall.entry[names(markers)]

    stat.order <- order(stats,decreasing=T)

    unique.pvals <- unique.pvals[stat.order]
    markers <- markers[stat.order]
    fcs <- fcs[stat.order] / log(2)
    stats <- stats[stat.order]
    total.fcs <- total.fcs[stat.order]

    marker.entry <- data.frame(sapply(markers,function(x) sprintf('%.2e',x)),
                               sapply(fcs,function(x) round(x,digits=2)),
                               sapply(stats,function(x) round(x,digits=2)),
                               sapply(total.fcs,function(x) round(x,digits=2)),
                               sapply(unique.pvals,function(x) sprintf('%.2e',x)))

    rownames(marker.entry) <- names(markers)
    marker.list[[i]] <- marker.entry
}

saveRDS(marker.list,'results/recount/markers/pairwise_leiden_tpmlogclus_pc3sd_poscounts_marker_list_filtered_8thpval_8thfc.RDS')
