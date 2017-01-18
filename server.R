library(plyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gtools)
library(stats)
library(shiny)
library(scales)
library(plotly)

shinyServer(function(input,output){  

    ## output$tsne <- renderPlot({
    output$tsne <- renderPlotly({
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

        tsne.y <- readRDS('data/plotData.RDS')
        
        
        ## tsne.y$kmeans.cluster <- as.factor(tsne.y$kmeans.cluster)

        tsne <- plot_ly(tsne.y,x=~y1,y=~y2,type='scatter',mode='markers',
                        hoverinfo = 'text',
                        text = ~paste(tissue_general,'<br>',project,'<br>',ethnicity))
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

        if(input$genevsgroup == 1){
            ## tsne <- tsne + geom_point(colour='black',size=2.9,alpha=0.9)
            ## tsne <- tsne %>%
            ##     add_trace(colour='black')
        }

        if(input$genevsgroup == 2){
            if(length(input$colorfactors) == 0){
                ## tsne <- tsne + geom_point(colour='black',size=2.9,alpha=0.9)
            }
            else{
                ## tsne <- tsne + geom_point(aes(colour=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+')),size=2.9,alpha=0.9) + guides(color=guide_legend(title='Metadata Group'))

                tsne <- tsne %>%
                add_trace(color=apply(tsne.y[,input$colorfactors,drop=FALSE],1,paste,collapse='+'))
            }
        }
        
        ## if(input$genevsgroup == 3){
        ##     if(length(input$groupShape) == 0){
        ##         tsne <- tsne + geom_point(aes(colour=tsne.y[,input$genecolor]),size=2.9,alpha=0.9) + scale_colour_gradient2(low ='red',mid='white',high='blue',limits=c(min(tsne.y[,input$genecolor]),max(tsne.y[,input$genecolor]))) + guides(color=guide_legend(title = sprintf('%s Expression',input$genecolor)))
        ##     }
        ##     else{
        ##         tsne <- tsne + geom_point(aes(colour=tsne.y[,input$genecolor],shape=tsne.y[,input$groupShape]),size=2.9,alpha=0.9) + scale_colour_gradient2(low = 'red',mid='white',high='blue',limits=c(min(tsne.y[,input$genecolor]),max(tsne.y[,input$genecolor]))) + guides(color=guide_legend(title = sprintf('%s Expression',input$genecolor)),shape=guide_legend(title=input$groupShape))
        ##     }
        ## }
        tsne
    })
})


