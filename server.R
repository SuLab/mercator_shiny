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
library(jsonlite)
library(DT)

## lasso2d = '{
##     name: 'lasso2d',
##     title: 'Lasso Select',
##     attr: 'dragmode',
##     val: 'lasso',
##     icon: Icons.lasso,
##     click: handleCartesian
## };'

shinyServer(function(input,output,session){  

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
    tsne.order <- rownames(tsne.data[[1]])

    tsne.meta <- readRDS('data/gtex_tcga_meta.RDS')
    rownames(tsne.meta) <- tsne.meta$data_id
    tsne.meta <- tsne.meta[tsne.order,]

    gene.tpm.samps <- readRDS('data/tpm_mat_rownames.RDS')

    ## gene.choices <- readRDS('data/gene_ids.RDS')
    gene.choices <- readRDS('data/gene_label_vector.RDS')

    kmeans.dat <- readRDS('data/kmeans_1to100_recount_250_dim_noScaled_noProj_over50_noSingle_entrez.RDS')

    ## louvain.vec <- readRDS('data/louvain_pca_over50_k30.RDS')
    ## louvain.vec <- readRDS('data/louvain_pca_over50_k100.RDS')
    ## louvain.vec <- readRDS('data/tcga_sklearn_louvain.RDS')
    louvain.vec <- readRDS('data/kmeans_2_louvain_recount_k5_over50.RDS')
    louvain.choices <- sort(unique(readRDS('data/kmeans_2_louvain_recount_k5_over50.RDS')))
    louvain.choices <- sapply(louvain.choices,function(x) sprintf('Louvain Cluster %s',x))
    names(louvain.choices) <- louvain.choices
    
    marker.list <- readRDS('data/pairwise_kmeans_k30_marker_list_filtered_noDropouts.RDS')
    ## marker.list <- readRDS('data/marker_test.RDS')

    user.selections <- reactiveValues(
        selection.list = list(),
        selection.datalist = data.frame(),
        gene.input.results = NULL,
        tsne.traces = list(tsne.order)
    )
    

    ## selection.list <- list()
    ## selection.datalist <- data.frame()

    updateSelectizeInput(session,
                         'whichGene',
                         choices=gene.choices,
                         server=TRUE,
                         options=list(maxOptions=10,
                                      closeAfterSelect=TRUE
                                      )
                         )


    observeEvent(input$kmeansChoice, {

        kmeans.num <- as.integer(input$kmeansChoice)

        updateSelectInput(session,
                          'barPlotFactor',
                          label='Bar plot group',
                          choices=c('All','Selections',sapply(1:kmeans.num,function(x) sprintf('Kmeans Cluster %d',x)),louvain.choices)
                          )
    })

    selection.getDatalist <- reactive({

        input$saveSelection

        ## print('get datalist')
        ## print(user.selections$selection.list)
        
        if(length(user.selections$selection.list) == 0){

            return(data.frame('name'='demo','length'=0))

        }
        else{
            return(user.selections$selection.datalist)
        }
    })

    observeEvent(input$saveSelection, {
        
        event.data <- event_data('plotly_selected',source = 'tsne')

        if(is.null(event.data) == TRUE) return(NULL)

        curveNumber <- event.data$curveNumber+1

        inds <- event.data$pointNumber+1

        user.selections$selection.list[[length(user.selections$selection.list)+1]] <- user.selections$tsne.traces[[curveNumber[1]]][inds]

        ## user.selections$selection.list[[length(user.selections$selection.list)+1]] <- data.frame(curveNumber,inds)

        ## print(length(user.selections$selection.list))

        label <- input$selectionName

        if(label == ''){
            label <- sprintf('Selection %d',length(user.selections$selection.list))
        }

        user.selections$selection.datalist <- rbind(user.selections$selection.datalist,data.frame('name'=label,'length'=nrow(event.data)))

    })


    output$markerTable <- DT::renderDT({

        dataDT <- DT::datatable(marker.list[[as.double(input$markerGroup)]],
                                colnames=c('P-val','log fc'),
                                class='display nowrap',
                                options = list(
                                    ## searching=TRUE,
                                    'lengthChange'=FALSE,
                                    'info'=FALSE,
                                    'pagingType'='simple')
                                )

        return(dataDT)
    })


    output$selectionList <- DT::renderDT({

        datalist <- selection.getDatalist()


        dataDT <- DT::datatable(datalist,
                                options = list(
                                    searching = FALSE,
                                    paging = FALSE,
                                    info = FALSE)
                                )
        
        ## print(datalist)
        ## print(selection.list)

        ## return(datalist)
        return(dataDT)
    })

    observeEvent(input$gene.vec, {

        req <- POST('http://localhost:3000/euclid_pca',body=upload_file('euclid.pca.file$datapath','text/plain'))
        stop_for_status(rep)

        req.text <- content(req,'text','text/plain')

        con <- textConnection(req.text)
        dist.vec <- read.table(con,sep='\t',header=T)

        print(head(dist.vec))

        dat.rows <- rownames(data())
        ## dat.rows <- gsub('-','.',dat.rows)
        ## dat.rows

        user.selections$gene.input.results <- dist.vec[dat.rows,]

    })


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


    geneVec <- reactive({

        ## geneVec <- NULL

        if(input$whichGene == ''){
            return(NULL)
        }
        
        gene.id <- strsplit(input$whichGene,', ')[[1]][3]

        gene.info <- fromJSON(sprintf('http://localhost:3000/gene_vals/%s',gene.id))

        names(gene.info) <- gene.tpm.samps

        geneVec <- rep(NA,length(tsne.order))
        names(geneVec) <- tsne.order

        names.used <- intersect(tsne.order,gene.tpm.samps)

        geneVec[names.used] <- gene.info[names.used]

        ## gene.info <- gene.info[tsne.order]

        return(geneVec)

    })
        

    meshVec <- reactive({

        tree <- input$tree

        if(is.null(tree) | length(get_selected(tree)) == 0){
            colVar=NULL
        } else{
            mesh.selection <- unlist(get_selected(tree))
            mesh.id <- strsplit(mesh.selection,':')[[1]][1]
            mesh.info <- fromJSON(sprintf('http://localhost:3000/ontology_info/%s',mesh.id))
            colVar <- rep('unlabelled',length(tsne.order))
            names(colVar) <- tsne.order

            for(label in names(mesh.info)){
                used.ids <- intersect(tsne.order,mesh.info[[label]])
                colVar[used.ids] <- label
            }
        }
        return(colVar)
    })

    selectionVec <- reactive({

        rows <- input$selectionList_rows_selected

        if(is.null(rows)){
            return(NULL)
        }

        indVec <- c()
        groupVec <- c()

        dataResults <- data()

        results <- data.frame()

        for(rowNum in rows){

            entry <- user.selections$selection.list[[rowNum]]

            rowLabel <- sprintf('Selection: %s',user.selections$selection.datalist[rowNum,'name'])

            ## results <- rbind(results,data.frame(samps = user.selections$tsne.traces[[entry$curveNumber[1]]][entry$inds],group=rep(rowLabel,nrow(entry))))

            results <- rbind(data.frame(samps=entry,group=rep(rowLabel,length(entry))))

            ## results <- rbind(results,data.frame(samps=rownames(dataResults)[entry],group=rep(rowLabel,length(entry))))

            ## indVec <- c(indVec,entry)
            ## groupVec <- c(groupVec,rep(rowNum,length(entry)))

        }

        

        ## return(list(
        ##     indices = rownames(dataResults)[indVec],
        ##     groups = groupVec
        ## ))

        return(results)

    })

    selectionBarVec <- reactive({

        rows <- input$selectionList_rows_selected

        if(is.null(rows)){
            rows <- 1:length(user.selections$selection.list)
        }

        indVec <- c()
        groupVec <- c()

        dataResults <- data()

        results <- data.frame()

        for(rowNum in rows){

            entry <- user.selections$selection.list[[rowNum]]

            rowLabel <- sprintf('Selection: %s',user.selections$selection.datalist[rowNum,'name'])

            ## results <- rbind(results,data.frame(samps = user.selections$tsne.traces[[entry$curveNumber[1]]][entry$inds],group=rep(rowLabel,nrow(entry))))

            results <- rbind(results,data.frame(samps=entry,group=rep(rowLabel,length(entry))))

            ## results <- rbind(results,data.frame(samps=rownames(dataResults)[entry],group=rep(rowLabel,length(entry))))

            ## indVec <- c(indVec,entry)
            ## groupVec <- c(groupVec,rep(rowNum,length(entry)))

        }

        

        ## return(list(
        ##     indices = rownames(dataResults)[indVec],
        ##     groups = groupVec
        ## ))

        return(results)

    })


    violinColVar <- reactive({
        colVar <- NULL
        if(input$colorButton == 0 | input$violinXFactors %in% c('No Coloring','Gene','Mesh','KMeans','Louvain')){
            return(colVar)
        }

        colVar = apply(data()[,input$violinXFactors,drop=FALSE],1,paste,collapse='+')

        return(colVar)

    })

    kmeansVec <- reactive({

        ## print(input$kmeansChoice)
        
        return(kmeans.dat[[as.numeric(input$kmeansChoice)]])

    })

    louvainVec <- reactive({

        return(louvain.vec)
    })

    colVar <- reactive({
        
        colVar <- NULL
        ## if(TRUE){return(colVar)}
        if(input$colorButton == 0 | input$colorfactors %in% c('No Coloring','Gene','Mesh','KMeans','Louvain')){
            return(colVar)
        }

        ## geneVecResults <- geneVec()

        ## if(length(input$colorfactors) > 0){

        ##     print(input$colorfactors)

        ## if(input$colorfactors == 'Mesh'){
            
        ##     tree <- input$tree

        ##     if(is.null(tree)){
        ##         colVar=NULL}
        ##     else{
        ##         mesh.selection <- unlist(get_selected(tree))
        ##         mesh.id <- strsplit(mesh.selection,':')[[1]][1]
        ##         mesh.info <- fromJSON(sprintf('http://localhost:3000/ontology_info/%s',mesh.id))
        ##         colVar <- rep('unlabelled',length(tsne.order))
        ##         names(colVar) <- tsne.order

        ##         for(label in names(mesh.info)){
        ##             used.ids <- intersect(tsne.order,mesh.info[[label]])
        ##             colVar[used.ids] <- label
        ##         }
        ##     }

        ## } else
        ## if(input$colorfactors == 'Gene') {

        ##     colVar <- geneVecResults
            
        ## } else{

        colVar = apply(data()[,input$colorfactors,drop=FALSE],1,paste,collapse='+')
            
        ## }

        return(colVar)
        ## }
    })

    barColVar <- reactive({
        
        colVar <- NULL
        ## if(TRUE){return(colVar)}
        if(input$colorButton == 0 | input$barPlotXaxis %in% c('No Coloring','Gene','Mesh','KMeans','Louvain')){
            return(colVar)
        }

        colVar = apply(data()[,input$barPlotXaxis,drop=FALSE],1,paste,collapse='+')
            
        return(colVar)

    })

    tree.dat <- readRDS('data/mesh_tree_flat.RDS')

    output$tree <- renderTree(tree.dat)

    output$metadataBar <- renderPlot({
        
        input$colorButton

        isolate({
            ## euclid.file <- input$euclid_input
            ## spear.file <- input$spearman_input
            ## euclid.pca.file <- input$euclid_pca_input
            colVarResults <- barColVar()
            dataResults <- data()
            meshResults <- meshVec()
            barCategory <- input$barPlotXaxis
            barGroup <- input$barPlotFactor
            selectionRes <- selectionBarVec()
            kmeansVec <- kmeansVec()
            louvainVec <- louvainVec()            
        })

        if(input$colorButton == 0 | barCategory == 'No Coloring'){
            xGroup <- 'All'
        } else if(barCategory == 'Mesh'){
            xGroup <- meshResults
        } else{
            xGroup <- colVarResults
        }

        saveRDS(selectionRes,'../../data/tmp/selection.RDS')

        if(barGroup == 'Selections'){
            xGroup <- xGroup[as.character(selectionRes$samps)]

        }else if(grepl('Kmeans',barGroup)){
            kmeans.group <- as.integer(gsub('Kmeans Cluster ','',barGroup))
            used.samps <- names(which(kmeansVec==kmeans.group))
            xGroup <- xGroup[as.character(used.samps)]

        }else if(grepl('Louvain',barGroup)){
            print('allo')
            louvain.group <- gsub('Louvain Cluster ','',barGroup)
            used.samps <- names(which(louvainVec==louvain.group))
            xGroup <- xGroup[as.character(used.samps)]
        }

        labels <- unlist(lapply(xGroup,function(x) strsplit(x,' [+] ')[[1]]))
        ## for(label in xGroup){
        ##     label.split <- strsplit(label,' [+] ')[[1]]
        ##     labels <- c(labels,label.split)
        ## }
        ## saveRDS(labels,'foo')
        ## for(label in unique(xGroup)){
            
        ##     label.split <- strsplit(label,' + ')[[1]]
        ##     labels <- c(labels,label.split)

        ## }

        metadataTable <- table(labels)
        ## metadataTable <- table(xGroup)

        plot.dat <- data.frame(Label = names(metadataTable), number = as.vector(metadataTable))
        
        if(xGroup == 'All'){
            
            plot.dat <- data.frame(Label= c('All'), number = c(length(tsne.order)))
        }
        
        output.plot <- ggplot() +
            geom_bar(data=plot.dat,aes(x=Label,y=number,fill=Label),stat='identity')

        output.plot

    })
        

    output$violin <- renderPlot({

        input$colorButton

        isolate({
            euclid.file <- input$euclid_input
            spear.file <- input$spearman_input
            euclid.pca.file <- input$euclid_pca_input
            colVarResults <- violinColVar()
            dataResults <- data()
            geneVecResults <- geneVec()
            meshResults <- meshVec()
            colorFactors <- input$violinXFactors
            kmeansVec <- kmeansVec()
            louvainVec <- louvainVec()
            selectionRes <- selectionVec()
        })

        ## plot.dat <- cbind(dataResults,colVarResults)
        ## colnames(plot.dat) <- c('y1','y2','color')

        if(input$colorButton == 0 | colorFactors == 'No Coloring'){
            xGroup <- 'All'
        } else if(colorFactors == 'Mesh'){
            xGroup <- meshResults
        } else if(colorFactors == 'KMeans'){
            xGroup <- kmeansVec
            ## xGroup <- colVarResults
            ## xGroup <- colorFactors
        } else if(colorFactors == 'Louvain'){
            xGroup <- louvainVec
        }
        else{
            xGroup <- colVarResults
        }
        ## if(is.null(xGroup)){
        ##     yGroup <- NULL
        ## }
        ## else{
        yGroup <- geneVecResults
        ## }

        ## print(head(names(yGroup)))
        ## print(head(names(xGroup)))

        ## saveRDS(selectionRes,'../../data/tmp/selection.RDS')

        if(any(xGroup != 'All')){
            used.names <- intersect(names(xGroup),names(yGroup))
            plot.dat <- data.frame(x=as.factor(xGroup[used.names]),y=yGroup[used.names])

            if(!is.null(selectionRes)){

                selectionRes <- selectionRes[selectionRes$samps %in% used.names,]

                plot.dat <- rbind(plot.dat,data.frame(x=as.factor(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))

            }
        }
        else{
            ## used.names <- names(yGroup)
            plot.dat <- data.frame(x=as.factor(xGroup),y=yGroup)

        }



        plot.dat <- subset(plot.dat,x != 'unlabelled')
        ## print(table(xGroup))
        output.plot <- ggplot() +
            geom_violin(data=plot.dat,aes(x=x,y=log(y+1),fill=x),scale='width') +
            geom_jitter(data=plot.dat,aes(x=x,y=log(y+1)),width=0.1,alpha=0.1,colour='gray') +
            theme(panel.background= element_blank(),
                  axis.text.x = element_text(size=24,angle=75,vjust=0.5),
                  axis.text.y = element_text(size=24),
                  legend.text = element_text(size=24),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size=28),
                  legend.title = element_blank(),
                  legend.position='none'
                  )
            ## geom_violin(data=dataResults, x = ~y1, y = ~y2,mode="markers",type='scatter',color = colVarResults,text=colVarResults)
            ## config(p = .,modeBarButtonsToRemove = c("zoom2d",'toImage','autoScale2d','hoverClosestGl2d'),collaborate=FALSE,cloud=FALSE) %>%
            ## config(collaborate=FALSE) %>%

        output.plot
        
    })

    output$tsne <- renderPlotly({

        input$colorButton

        isolate({
            euclid.file <- input$euclid_input
            spear.file <- input$spearman_input
            euclid.pca.file <- input$euclid_pca_input
            colVarResults <- colVar()
            dataResults <- data()
            geneVecResults <- geneVec()
            meshResults <- meshVec()
            colorFactors <- input$colorfactors
            kmeansVec <- kmeansVec()
            louvainVec <- louvainVec()
            ## input$tree
            ## input$colorfactors
        })


        ## print(colorFactors)
        
        if(input$colorButton == 0 | colorFactors == 'No Coloring'){
            colVarPlot <- NULL

        } else if(colorFactors == 'Mesh'){
            colVarPlot <- meshResults

        } else if(colorFactors == 'Gene'){
            colVarPlot <- geneVecResults

        } else if(colorFactors == 'KMeans'){
            colVarPlot <- as.factor(kmeansVec[rownames(dataResults)])

        } else if(colorFactors == 'Louvain'){
            colVarPlot <- as.factor(louvainVec[rownames(dataResults)])

        }
        else{
            colVarPlot <- colVarResults
        }


        if(!is.null(colVarPlot) & colorFactors != 'Gene'){
            
            i <- 1

            ## isolate({user.selections$tsne.traces <- list()})

            tsne.traces <- list()

            for(label in sort(unique(colVarPlot))){

                inds <- which(colVarPlot == label)

                tsne.traces[[i]] <- rownames(dataResults)[inds]
                
                ## print(head(user.selections$tsne.traces[[i]]

                i <- i+1
            }

            isolate({
                user.selections$tsne.traces <- tsne.traces
            })

        }
        else{
            isolate({
                user.selections$tsne.traces <- list(rownames(dataResults))
            })
        }

        ## print(colVarPlot)
        ## print(typeof(colVarPlot))

        ## print('hello')

        output.plot <- plot_ly(dataResults, x = ~y1, y = ~y2,mode="markers",type='scatter',color = colVarPlot,
                               text=colVarPlot,
                               source='tsne') %>%
            ## config(p = .,modeBarButtonsToRemove = c("zoom2d",'toImage','autoScale2d','hoverClosestGl2d'),collaborate=FALSE,cloud=FALSE) %>%
            ## config(collaborate=FALSE) %>%
            layout(dragmode = "pan",xaxis=ax,yaxis=ax) %>%
            toWebGL()

        output.plot
        

        ## output.plot <- plot_ly(data(), x = ~y1, y = ~y2,mode="markers",type='scatter',color = colVar,hoverInfo='text',text=colVar) %>%
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






### stuff for euclid

        ## colVar <- NULL
        ## ## colVar = NULL
        ## ## if(input$genevsgroup == 0){
        ## ## }

        ## ## if(input$genevsgroup == 2){
        ## ## if(length(input$colorfactors) == 0){
            
        ## ## }
        ## if(input$colorFactors == 'Mesh'){
            
        ##     tree <- input$tree

        ##     if(is.null(tree)){
        ##         colVar=NULL}
        ##     else{
        ##         mesh.selection <- unlist(get_selected(tree))
        ##         mesh.id <- strsplit(mesh.selection,':')[[1]][1]
        ##         mesh.info <- fromJSON(sprintf('http://localhost:3000/ontology_info/%s',mesh.id))
        ##         colVar <- rep('unlabelled',length(tsne.order))
        ##         names(colVar <- tsne.order)
        ##         for(label in names(mesh.info)){
        ##             colVar[mesh.info[[label]]] <- label
        ##         }
        ##     }

        ## } else if(length(input$colorfactors) > 0){

        ##     colVar = apply(data()[,input$colorfactors,drop=FALSE],1,paste,collapse='+')
            
        ## }

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
