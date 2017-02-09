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

shinyServer(function(input,output){  

    ## data <- reactive({
    ##     tsne.y <- readRDS('data/plotData.RDS')
    ##     return(tsne.y)
    ## })

    data <- reactive({
        tsne.y <- readRDS('data/plotData.RDS')
        ## rownames(tsne.y) <- NULL

        ## mod.y <- tsne.y

        ## for(i in 1:4){
        ##     mod.y$y1 <- runif(nrow(mod.y),-60,60)
        ##     mod.y$y2 <- runif(nrow(mod.y),-60,60)

        ##     tsne.y <- rbind(tsne.y,mod.y)
        ## }
        
        
        return(tsne.y)
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
    output$tsne <- renderPlotly({
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
        colVar = NULL
        if(input$genevsgroup == 1){
            ## tsne <- tsne + geom_point(colour='black',size=2.9,alpha=0.9)
            ## tsne <- tsne %>%
            ##     add_trace(colour='black')
        }


        ## if(!(is.null(input$gene.vec))){
        ##     input.dat <- read.table(input$gene.vec$datapath,sep='\t',header=F,stringsAsFactors=F)

        ##     head(input.dat)
        ##     dist.vec <- apply(tsne.y[,rownames(input.dat)],1,dist)

        ##     tsne.y$distances <- dist.vec
            
        ##     tsne <- tsne %>%
        ##         add_trace(color=~distances)
        ## }
        ## else{
        if(input$genevsgroup == 2){
            if(length(input$colorfactors) == 0){
                ## tsne <- tsne + geom_point(colour='black',size=2.9,alpha=0.9)
            }
            else{
                ## tsne <- tsne + geom_point(aes(colour=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+')),size=2.9,alpha=0.9) + guides(color=guide_legend(title='Metadata Group'))

                ## tsne <- tsne %>%
                ##     add_trace(color=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+'))
                ## tsne <- scatterD3(data=tsne.y,x=y1,y=y2,col_var=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+'),transitions=TRUE)

                colVar = apply(data()[,input$colorfactors,drop=FALSE],1,paste,collapse='+')
                
            }
        }
        ## }
        
        ## if(input$genevsgroup == 3){
        ##     if(length(input$groupShape) == 0){
        ##         tsne <- tsne + geom_point(aes(colour=tsne.y[,input$genecolor]),size=2.9,alpha=0.9) + scale_colour_gradient2(low ='red',mid='white',high='blue',limits=c(min(tsne.y[,input$genecolor]),max(tsne.y[,input$genecolor]))) + guides(color=guide_legend(title = sprintf('%s Expression',input$genecolor)))
        ##     }
        ##     else{
        ##         tsne <- tsne + geom_point(aes(colour=tsne.y[,input$genecolor],shape=tsne.y[,input$groupShape]),size=2.9,alpha=0.9) + scale_colour_gradient2(low = 'red',mid='white',high='blue',limits=c(min(tsne.y[,input$genecolor]),max(tsne.y[,input$genecolor]))) + guides(color=guide_legend(title = sprintf('%s Expression',input$genecolor)),shape=guide_legend(title=input$groupShape))
        ##     }
        ## }

        output.plot <- plot_ly(data(), x = ~y1, y = ~y2) %>% toWebGL()
        output.plot
        ## scatterD3(x = data()[,'y1'],
        ##           y = data()[,'y2'],
        ##           col_var = colVar)
                  ## transitions=TRUE)
    })
})
