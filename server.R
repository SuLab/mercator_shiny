library(plyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gtools)
library(stats)
library(shiny)
library(scales)
library(plotly)
library(scatterD3)
library(httr)
library(shinyTree)


## lasso2d = '{
##     name: 'lasso2d',
##     title: 'Lasso Select',
##     attr: 'dragmode',
##     val: 'lasso',
##     icon: Icons.lasso,
##     click: handleCartesian
## };'

shinyServer(function(input,output){  

    ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    )

    ## tsne.data <- readRDS('data/recount_tsne_pca_list_noNAs_tcgagtex.RDS')
    ## tsne.data <- readRDS('data/recount_tsne_pca_list_noNAs_tcgagtex_fewer_genes.RDS')
    ## tsne.data <- readRDS('data/recount_tsne_pca_list_noNAs_tcgagtex.RDS')
    tsne.data <- readRDS('data/recount_tsne_pca_list_noNAs_over50_noSingle_entrez.RDS')


    tsne.meta <- readRDS('data/gtex_tcga_meta.RDS')
    rownames(tsne.meta) <- tsne.meta$data_id
    tsne.meta <- tsne.meta[rownames(tsne.data[[1]]),]

    data <- reactive({
        ## tsne.data <- readRDS('data/recount_tsne.RDS')
        ## tsne.data <- readRDS('data/recount_tsne_pca_list_noNAs_tcgagtex.RDS')
        ## tsne.y <- data.frame(y1=tsne.data$Y[,1],y2=tsne.data$Y[,2])

        tsne.y <- data.frame(tsne.data[[input$whichTSNE]])
        colnames(tsne.y) <- c('y1','y2')

        tsne.y <- cbind(tsne.y,tsne.meta)

        ## tsne.meta <- as.data.frame(readRDS('data/all_recount_metasra_summarized.RDS'))
        ## tsne.meta <- tsne.meta[rownames(tsne.y),]

        ## tsne.y <- cbind(tsne.y,tsne.meta)
        ## creates extra datapoints
        ## rownames(tsne.y) <- NULL

        ## mod.y <- tsne.y

        ## for(i in 1:4){
        ##     mod.y$y1 <- runif(nrow(mod.y),-60,60)
        ##     mod.y$y2 <- runif(nrow(mod.y),-60,60)

        ##     tsne.y <- rbind(tsne.y,mod.y)
        ## }
        
        ## tsne.y$tissue_detail[unlist(sapply(tsne.y$tissue_detail,function(x) !grepl("Adenocarcinoma",x)))] <- "NA"
        ## tsne.y$tissue_general[tsne.y$tissue_general != "hESC"] <- "NA"

        ## return(tsne.y)
        return(tsne.y)
    })

    tree.dat <- readRDS('data/mesh_tree_flat.RDS')

    output$tree <- renderTree(tree.dat)
        

    output$tsne <- renderPlotly({

        euclid.file <- input$euclid_input
        spear.file <- input$spearman_input
        euclid.pca.file <- input$euclid_pca_input

        colVar <- NULL
        ## colVar = NULL
        ## if(input$genevsgroup == 0){
        ## }

        ## if(input$genevsgroup == 2){
        ## if(length(input$colorfactors) == 0){
            
        ## }
        if(length(input$colorfactors) > 0){

            colVar = apply(data()[,input$colorfactors,drop=FALSE],1,paste,collapse='+')
            
        }
        ## }

        ## if(!is.null(spear.file)){

        ##     req <- POST("http://localhost:3000/spearman_distance",body=upload_file(spear.file$datapath,'text/plain'))
        ##     stop_for_status(req)
            
        ##     req.text <- content(req,"text","text/plain")

        ##     con <- textConnection(req.text)
        ##     dist.vec <- read.table(con,sep='\t',header=T)

        ##     dat.rows <- rownames(data())
        ##     dat.rows <- gsub('-','.',dat.rows)
        ##     dat.rows <- sub('^([0-9])','X\\1',dat.rows)

        ##     colVar <- dist.vec[dat.rows,]

        ##     colVar <- order(colVar)

        ##     ## print(min(colVar))

        ##     ## ## colVar[colVar==min(colVar)] <- min(subset(dist.vec,x>1))

        ##     ## print(min(colVar))
        ## }

        ## if(!is.null(euclid.file)){

        ##     req <- POST("http://localhost:3000/euclid_distance",body=upload_file(euclid.file$datapath,'text/plain'))
        ##     stop_for_status(req)
            
        ##     req.text <- content(req,"text","text/plain")

        ##     con <- textConnection(req.text)
        ##     dist.vec <- read.table(con,sep='\t',header=T)

        ##     dat.rows <- rownames(data())
        ##     dat.rows <- gsub('-','.',dat.rows)
        ##     dat.rows <- sub('^([0-9])','X\\1',dat.rows)

        ##     colVar <- dist.vec[dat.rows,]

        ##     colVar[colVar==min(colVar)] <- min(subset(dist.vec,x>1))

        ##     colVar <- order(colVar)

        ## }

        ## if(!is.null(euclid.pca.file)){

        ##     req <- POST("http://localhost:3000/euclid_pca",body=upload_file(euclid.pca.file$datapath,'text/plain'))
        ##     stop_for_status(req)
            
        ##     req.text <- content(req,"text","text/plain")

        ##     con <- textConnection(req.text)
        ##     dist.vec <- read.table(con,sep='\t',header=T)

        ##     print(head(dist.vec))

        ##     dat.rows <- rownames(data())
        ##     dat.rows <- gsub('-','.',dat.rows)
        ##     dat.rows <- sub('^([0-9])','X\\1',dat.rows)

        ##     colVar <- dist.vec[dat.rows,]
        ##     ## colVar <- sqrt(colVar)
        ##     ## colVar <- max(colVar) - colVar


        ##     colVar[colVar==min(colVar)] <- min(subset(dist.vec,x>1))
        ##     colVar <- log(colVar)
        ##     colVar <- 1.01 ^ (max(colVar) - colVar)

        ##     ## colVar <- order(colVar)

        ## }
        

        ## output.plot <- plot_ly(data(), x = ~y1, y = ~y2,mode="markers",type='scatter',color = colVar,hoverInfo='text',text=colVar) %>%
        output.plot <- plot_ly(data(), x = ~y1, y = ~y2,mode="markers",type='scatter',color = colVar,text=colVar) %>%
            ## config(p = .,modeBarButtonsToRemove = c("zoom2d",'toImage','autoScale2d','hoverClosestGl2d'),collaborate=FALSE,cloud=FALSE) %>%
            ## config(collaborate=FALSE) %>%
            layout(dragmode = "pan",xaxis=ax,yaxis=ax) %>%
            toWebGL()

        output.plot
    })
})

    ## tsne.y <- readRDS('data/foo.RDS')

    ## tsne.y <- readRDS('data/plotData.RDS')
    
    ## gene.dat <- read.table('data/gene_dat.tsv',sep='\t',header=T,stringsAsFactors=F)

    ## tsne.y <- cbind(tsne.y,gene.dat)

    ## tsne.y$kmeans.cluster <- as.factor(tsne.y$kmeans.cluster)

    ## plotDat <- reactive({
    ##     readRDS('data/plotData.RDS')
    ##     })


    ## output$tsne <- renderPlot({



    ## output$tsne <- renderScatterD3({

        ## tsne.y <- plotDat()
        
        ## ctDR <- read.table('data/scqpcr.csv',sep=',',header=T)

        ## ## remove duplicated rows
        ## if(anyDuplicated(ctDR[,10:ncol(ctDR)]) > 0){
        ##     ctDR <- ctDR[-which(duplicated(ctDR[,10:ncol(ctDR)])),]
        ## }

        ## ## remove genes without variance
        ## vars <- NULL
        ## for(i in 10:ncol(ctDR)){
        ##     vars <- c(vars, var(ctDR[,i], na.rm=T))
        ## }
        ## ctGenesNoVar <- ctDR[which(vars == 0 | is.na(vars))+9]
        ## ctDR <- ctDR[,c(1:9, (which(!is.na(vars) & vars!=0)+9))]


        ## ## get ordered gene list from differential expression between clusters
        ## ## get p-values for differential expression for factors
        ## pvals <- NULL
        ## for(i in 10:ncol(ctDR)){
        ##     pvals <- c(pvals, kruskal.test(ctDR[,i], factor(ctDR[, "kmeans.cluster"]))$p.value)
        ## }

        ## ctDR <- ctDR[,c(1:9, (order(pvals)+9))]
        ## genesByDiff <- names(ctDR)[10:ncol(ctDR)]


        ## ## t-SNE
        ## set.seed(1)
        ## tsne_out <- Rtsne(as.matrix(ctDR[,10:ncol(ctDR)]), perplexity = 30)

        ## tsne_y <- as.data.frame(cbind(tsne_out$Y, 
        ##                               ctDR$cellSource,
        ##                               ctDR$kmeans.cluster,
        ##                               ctDR[, genesByDiff]))

        ## names(tsne_y)[1:4] <- c("y1", "y2", "tissue", "kmeans.cluster")
        ## kmeans.levels <- unique(tsne_y$kmeans.cluster)
        ## kmeans.levels <- kmeans.levels[order(kmeans.levels)]
        ## tsne_y$kmeans.cluster <- factor(tsne_y$kmeans.cluster, levels=kmeans.levels)

        ## for(i in c(1,2,5:ncol(tsne_y))){
        ##     tsne_y[, i] <- as.numeric(tsne_y[, i])
        ## }

        

        ## tsne <- plot_ly(tsne.y,x=~y1,y=~y2,type='scatter',mode='markers',
        ##                 hoverinfo = 'text',
        ##                 text = ~paste(tissue_general,'<br>',project,'<br>',tumor_sample_type))


        ## tsne <- scatterD3(data=tsne.y,x=y1,y=y2,transitions=TRUE)
        
        ## tsne <- ggplot(tsne.y, aes(y1,y2)) +
        ##     theme_bw() +
        ##     scale_x_continuous(breaks=seq(min(tsne.y$y1), max(tsne.y$y1), length.out = 10),
        ##                        minor_breaks = NULL) +
        ##     scale_y_continuous(breaks=seq(min(tsne.y$y2), max(tsne.y$y2), length.out = 10),
        ##                        minor_breaks = NULL) +
        ##     theme(axis.line=element_blank(),
        ##           axis.text.x=element_blank(),
        ##           axis.text.y=element_blank(),
        ##           axis.ticks=element_blank(),
        ##           axis.title.x=element_blank(),
        ##           axis.title.y=element_blank(),
        ##           panel.border=element_blank(),
        ##           panel.grid.major=element_blank(),
        ##           panel.grid.minor=element_blank()) 




            ## tsne <- tsne + geom_point(colour='black',size=2.9,alpha=0.9)
            ## tsne <- tsne %>%
            ##     add_trace(colour='black')

                ## tsne <- tsne + geom_point(colour='black',size=2.9,alpha=0.9)

    ## data <- reactive({
    ##     tsne.y <- readRDS('data/plotData.RDS')
    ##     return(tsne.y)
    ## })


        ## if(!(is.null(input$gene.vec))){
        ##     input.dat <- read.table(input$gene.vec$datapath,sep='\t',header=F,stringsAsFactors=F)

        ##     head(input.dat)
        ##     dist.vec <- apply(tsne.y[,rownames(input.dat)],1,dist)

        ##     tsne.y$distances <- dist.vec
            
        ##     tsne <- tsne %>%
        ##         add_trace(color=~distances)
        ## }
        ## else{


                ## tsne <- tsne + geom_point(aes(colour=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+')),size=2.9,alpha=0.9) + guides(color=guide_legend(title='Metadata Group'))

                ## tsne <- tsne %>%
                ##     add_trace(color=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+'))
                ## tsne <- scatterD3(data=tsne.y,x=y1,y=y2,col_var=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+'),transitions=TRUE)


        ## }
        
        ## if(input$genevsgroup == 3){
        ##     if(length(input$groupShape) == 0){
        ##         tsne <- tsne + geom_point(aes(colour=tsne.y[,input$genecolor]),size=2.9,alpha=0.9) + scale_colour_gradient2(low ='red',mid='white',high='blue',limits=c(min(tsne.y[,input$genecolor]),max(tsne.y[,input$genecolor]))) + guides(color=guide_legend(title = sprintf('%s Expression',input$genecolor)))
        ##     }
        ##     else{
        ##         tsne <- tsne + geom_point(aes(colour=tsne.y[,input$genecolor],shape=tsne.y[,input$groupShape]),size=2.9,alpha=0.9) + scale_colour_gradient2(low = 'red',mid='white',high='blue',limits=c(min(tsne.y[,input$genecolor]),max(tsne.y[,input$genecolor]))) + guides(color=guide_legend(title = sprintf('%s Expression',input$genecolor)),shape=guide_legend(title=input$groupShape))
        ##     }
        ## }




        ## scatterD3(x = data()[,'y1'],
        ##           y = data()[,'y2'],
        ##           col_var = colVar)
                  ## transitions=TRUE)
