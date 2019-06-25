library(plyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gtools)
library(stats)
library(shiny)
library(scales)
library(plotly)
library(httr)
library(shinyTree)
library(jsonlite)
library(DT)
library(shinyjs)
library(stringr)
library(parallel)
library(DBI)
library(odbc)

shinyServer(function(input,output,session){  

    ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        range = c(-55,55)
    )

    ## mercator.db.con <- dbConnect(odbc(),'PostgreSQL')

    tsne.data <- readRDS('data/recount_tsne_pca_noNAs_over50_bulkOnly_pc3sd_filteredp40.RDS')
    tsne.data <- tsne.data + matrix(data=rnorm(2*nrow(tsne.data),sd=0.5),ncol=2)    

    tsne.order <- rownames(tsne.data)

    tsne.meta <- readRDS('data/recount_meta_sampletype.RDS')
    rownames(tsne.meta) <- tsne.meta$run_id
    tsne.meta <- tsne.meta[tsne.order,]

    tcga.gtex.meta <- readRDS('data/gtex_tcga_meta.RDS')

    ## print(head(tcga.gtex.meta))
    
    rownames(tcga.gtex.meta) <- tcga.gtex.meta$data_id
    tcga.gtex.meta <- tcga.gtex.meta[tsne.order,]

    ## print(head(tcga.gtex.meta))

    tsne.meta <- cbind(tsne.meta,tcga.gtex.meta[,c('tissue_general','tissue_detail')])

    tcga.case.ids <- readRDS('data/pheno_case_ids.RDS')
    rownames(tcga.case.ids) <- tcga.case.ids$bigwig_file
    
    tcga.case.ids <- tcga.case.ids[tsne.order,]
    tsne.meta <- cbind(tsne.meta,tcga.case.ids$gdc_cases.case_id)

    gene.tpm.samps <- readRDS('data/tpm_mat_rownames.RDS')

    gene.choices <- readRDS('data/shiny_gene_choice_list.RDS')
    gene.num.vec <- readRDS('data/shiny_gene_num_vec.RDS')
    louvain.vec <- readRDS('data/leiden_pc3sd_r25e-3_vec.RDS')
    louvain.choices <- sort(unique(louvain.vec))
    louvain.choices <- sapply(louvain.choices,function(x) sprintf('Louvain Cluster %s',x))
    names(louvain.choices) <- louvain.choices

    ## top.clus.per.genes <- readRDS('data/leiden_pc3sd_r25e-3_markers_top_clus.RDS')


    marker.list <- readRDS('data/pairwise_leiden_r25e-3_marker_list_filtered_5th_perc_shiny.RDS')

    cellmarker.annotation.counts <- readRDS('data/cell_marker_annotations_marker_counts.RDS')

    user.selections <- reactiveValues(
        selection.list = list(),
        selection.datalist = data.frame(),
        gene.input.results = list(),
        gene.input.datalist = data.frame(),
        tsne.traces = list(tsne.order),
        previous.tab = 't-SNE'
    )


    updateSelectizeInput(session,
                         'whichGene',
                         choices=gene.choices,
                         server=TRUE,
                         options=list(maxOptions=10,
                                      closeAfterSelect=TRUE
                                      )
                         )



    projection.getDatalist <- reactive({

        input$gene.vec

        if(length(user.selections$gene.input.results) == 0){
            ## return(data.frame('name'=c(),'min'=c()))
            return(subset(data.frame('name'=c(1,2),'min'=c(2,3)),name>3))
        }
        else{
            return(user.selections$gene.input.datalist)
        }
    })

    selection.getDatalist <- reactive({

        input$saveSelection

        ## print('get datalist')
        ## print(user.selections$selection.list)
        
        if(length(user.selections$selection.list) == 0){

            ## return(data.frame('name'='demo','length'=0))
            ## return(data.frame('name'=c(),'length'=c()))
            return(subset(data.frame('name'=c(1,2),'length'=c(2,3)),name>3))

        }
        else{
            return(user.selections$selection.datalist)
        }
    })

    observeEvent(input$barPlotFactor, {

        if(grepl('Louvain', input$barPlotFactor)){
            updateSelectInput(session,
                              'barPlotXaxis',
                              label='X axis for Bar Plot',                              
                              choices = c('Project ID' = 'proj_id',
                                          'GTEx/TCGA Gross Tissue' = 'tissue_general',
                                          'GTEx/TCGA Detailed Tissue' = 'tissue_detail',
                                          'Sample Type' = 'sample_type',
                                          'Mesh: Anatomy' = 'Mesh',
                                          'Tissue' = 'Tissue',
                                          'DOID' = 'DOID',
                                          'EFO: Cultured Cells'='efo',
                                          'Louvain' = 'Louvain',
                                          'Marker Gene Tissue Types'='marker_tissue_type',
                                          'Marker Gene Cancer Type' = 'marker_cancer_type',
                                          'Marker Gene Cell Type' = 'marker_cell_type',
                                          'Marker Gene Cell Name' = 'marker_cell_name')
                              )
        } else{
            updateSelectInput(session,
                              'barPlotXaxis',
                              label='X axis for Bar Plot',                              
                              choices = c('No Coloring' = 'No Coloring',
                                          'Project ID' = 'proj_id',
                                          'GTEx/TCGA Gross Tissue' = 'tissue_general',
                                          'GTEx/TCGA Detailed Tissue' = 'tissue_detail',
                                          'Sample Type' = 'sample_type',
                                          'Mesh: Anatomy' = 'Mesh',
                                          'Tissue' = 'Tissue',
                                          'DOID' = 'DOID',
                                          'EFO: Cultured Cells'='efo',
                                          'Louvain' = 'Louvain')
                              )
        }
    })
                              
                              
        

    observeEvent(input$geneTable_rows_selected, {

        if(input$geneGroup == 'all'){
            selected.marker <- rownames(marker.list[['all']])[input$geneTable_rows_selected]
        }
        else{
            selected.marker <- rownames(marker.list[[as.double(input$geneGroup)+1]])[input$geneTable_rows_selected]
        }

        ## print(selected.marker)

        updateSelectizeInput(session,
                             'whichGene',
                             choices=gene.choices,
                             server=T,
                             options=list(maxOptions=10,
                                          closeAfterSelect=T
                                          ),
                             selected=gene.choices[gene.num.vec[selected.marker]]
                             )
        

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

    output$sampleInputTable <- DT::renderDT({

        datalist <- projection.getDatalist()

        dataDT <- DT::datatable(datalist,
                                selection='single',
                                options= list(
                                    searching=FALSE,
                                    paging=FALSE,
                                    info=FALSE)
                                )
        return(dataDT)
    })
                                                                
    output$geneTable <- DT::renderDT({


        if(input$geneGroup == 'all'){
            dt <- marker.list[['all']]
        } else{
            dt <- marker.list[[as.double(input$geneGroup)+1]]
        }

        dataDT <- DT::datatable(dt,
                                ## colnames=c('Symbol','ID','ENS_ID','gene-type','tissueType','cancerType','cellType','cellName','cluster p-val','cluster fc','cluster-all fc'),
                                selection='single',
                                options=list(
                                    'searching'=TRUE,
                                    'legnthChange'=FALSE,
                                    'pagingType'='full',
                                    'autoWidth'=TRUE,
                                    'pageLength'=25)
                                )
        

        return(dataDT)
    })

    output$sampleTable <- DT::renderDT({

        ## input$colorButton

        ## isolate({

        selectionRes <- selectionVec()
        louvainVec <- louvainVec()
        meshResults <- meshVec()
        tissueResults <- tissueVec()
        doidResults <- doidVec()
        efoResults <- efoVec()

        ## })

        dt <- tsne.meta[,c(1,2,3,4,5,6)]

        dt$Louvain <- louvainVec

        if(is.null(meshResults)){
            dt$mesh <- NA
        } else{
            dt$mesh <- meshResults
        }

        if(is.null(tissueResults)){
            dt$tissue <- NA
        } else{
            dt$tissue <- tissueResults
        }

        if(is.null(doidResults)){
            dt$doid <- NA
        } else{
            dt$doid <- doidResults
        }

        if(is.null(efoResults)){
            dt$efo <- NA
        } else{
            dt$efo <- efoResults
        }

        if(input$sampleTableControl == 'Selection'){
            dt <- dt[selectionRes$samps,]
            dt$group <- selectionRes$group
        ## } else if('Louvain' %in% input$sampleTableControl){
        } else if(grepl('Louvain',input$sampleTableControl)){
            louvain.id <- gsub('Louvain Cluster ','',input$sampleTableControl)

            dt <- subset(dt,Louvain==as.numeric(louvain.id))
        }

        rownames(dt) <- c()

        dataDT <- DT::datatable(dt,
                                options = list(pageLength=25
                                    )
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
        return(dataDT)
    })

    observeEvent(input$topNavBar, {

        if(input$topNavBar == 'genePanel'){
            
            runjs(" var $this = $('#controlPanelClickable');
                    if(!$this.hasClass('panel-collapsed')) {
                        $this.parents('.panel').find('.panel-body').slideUp();
		        $this.addClass('panel-collapsed');
		        $this.find('i').removeClass('glyphicon-chevron-up').addClass('glyphicon-chevron-down');
                    }")

        }
        else if(user.selections$previous.tab == 'genePanel'){

            runjs(" var $this = $('#controlPanelClickable');
                    if($this.hasClass('panel-collapsed')) {
                        $this.parents('.panel').find('.panel-body').slideDown();
		        $this.removeClass('panel-collapsed');
		        $this.find('i').removeClass('glyphicon-chevron-down').addClass('glyphicon-chevron-up');
                    }")
            
        }
        user.selections$previous.tab = input$topNavBar
    })

    observeEvent(input$gene.vec, {

        getDistance <- function(x,gene.vec) {

            eu.dist <- sum((x - gene.vec)^2)

            ## eu.dist <- (x %*% gene.vec) / (norm(x,'2') * norm(gene.vec,'2'))

            return(eu.dist)

        }


        progress <- shiny::Progress$new()

        on.exit(progress$close())

        progress$set(message='Calculating...', value=0.33)

        ## saveRDS(input$gene.vec,'../../data/tmp/input.RDS')

        map.dat <- readRDS('data/recount_750_dim_noProj_over50_bulkOnly_pc3sd_filtered.RDS')
        rotation.dat <- readRDS('data/recount_all_big_irlba_750_over50_bulkOnly_pc3sd_filtered_rotmat.RDS')
        map.centers <- readRDS('data/recount_over50_bulkOnly_pc3sd_filtered_colparams.RDS')
        map.centers <- subset(map.centers,vars>0)

        gene.dat <- read.table(input$gene.vec$datapath,sep='\t',header=T)
        gene.dat <- gene.dat[rownames(map.centers),,drop=F]

        dist.mat <- matrix(0,nrow=nrow(map.dat),ncol=ncol(gene.dat))
        rownames(dist.mat) <- rownames(map.dat)
        colnames(dist.mat) <- colnames(gene.dat)

        no_cores <- detectCores() - 1
        cl <- makeCluster(no_cores)

        for(col in colnames(gene.dat)){
    

            gene.vec <- gene.dat[[col]]
            names(gene.vec) <- rownames(gene.dat)

            center.vec <- gene.vec - map.centers$means
            center.vec <- center.vec / map.centers$vars

            rot.vec <- t(center.vec) %*% rotation.dat

            dist.vec <- apply(map.dat,1,getDistance,gene.vec=rot.vec)

            ## dist.vec <- parApply(cl,map.dat,1,getDistance,gene.vec=rot.vec)

            dist.mat[,col] <- dist.vec
        }

        stopCluster(cl)

        progress$set(message='Calculated', value=0.67)


        dat.rows <- rownames(data())

        dist.mat <- dist.mat[dat.rows,,drop=F]

        for(col in colnames(dist.mat)){

            ## dist.vec <- dist.mat[[col]]
            dist.vec <- dist.mat[,col]
            names(dist.vec) <- rownames(dist.mat)

            user.selections$gene.input.results[[length(user.selections$gene.input.results)+1]] <- dist.vec
            user.selections$gene.input.datalist <- rbind(user.selections$gene.input.datalist,data.frame('name'=col,'min'=min(log2(dist.vec+1))))
        }

        progress$set(message='Read results',value=1)

    })


    data <- reactive({
        tsne.y <- data.frame(tsne.data)
        colnames(tsne.y) <- c('y1','y2')

        tsne.y <- cbind(tsne.y,tsne.meta)

        return(tsne.y)
    })


    geneVec <- reactive({

        ## geneVec <- NULL

        ## if(input$whichGene == ''){
        if(TRUE){
            return(NULL)
        }
        
        ## gene.id <- strsplit(input$whichGene,', ')[[1]][3]

        gene.id <- input$whichGene

        gene.info <- fromJSON(sprintf('http://localhost:3000/gene_vals/%s',gene.id))

        names(gene.info) <- gene.tpm.samps

        geneVec <- rep(NA,length(tsne.order))
        names(geneVec) <- tsne.order

        names.used <- intersect(tsne.order,gene.tpm.samps)

        geneVec[names.used] <- gene.info[names.used]

        ## gene.info <- gene.info[tsne.order]

        ## return(geneVec)
        return(geneVec)
    })
        

    tissueVec <- reactive({

        tree <- input$tissueTree

        if(is.null(tree) | length(get_selected(tree)) == 0){
            colVar <- NULL
        } else{
            tissue.selection <- unlist(get_selected(tree))
            tissue.id <- strsplit(tissue.selection,':')[[1]][1]
            tissue.info <- fromJSON(sprintf('http://localhost:3000/tissue_info/%s',tissue.id))

            ## tissue.req <- dbFetch(dbSendQuery(mercator.db.con,sprintf("SELECT termtree FROM tissue_tree WHERE id='%s'",tissue.id)))
            ## tissue.info <- fromJSON(tissue.req$termtree)
            
            colVar <- rep('unlabelled',length(tsne.order))
            names(colVar) <- tsne.order

            for(label in names(tissue.info)){
                used.ids <- intersect(tsne.order,tissue.info[[label]])
                colVar[used.ids] <- label
            }
        }
        return(colVar)
    })

    doidVec <- reactive({

        tree <- input$doidTree

        if(is.null(tree) | length(get_selected(tree)) == 0){
            colVar <- NULL
        } else{
            doid.selection <- unlist(get_selected(tree))
            doid.split <- strsplit(doid.selection,':')[[1]]
            doid.id <- paste(doid.split[1],doid.split[2],sep='_')

            ## doid.req <- dbFetch(dbSendQuery(mercator.db.con,sprintf("SELECT termtree FROM doid_dat WHERE id='%s'",doid.id)))

            doid.info <- fromJSON(sprintf('http://localhost:3000/doid_info/%s',doid.id))
            ## doid.info <- fromJSON(doid.req$termtree)
            colVar <- rep('unlabelled',length(tsne.order))
            names(colVar) <- tsne.order

            for(label in names(doid.info)){
                used.ids <- intersect(tsne.order,doid.info[[label]])
                colVar[used.ids] <- label
            }
        }
        return(colVar)
    })

    efoVec <- reactive({

        tree <- input$efoTree

        if(is.null(tree) | length(get_selected(tree)) == 0){
            colVar <- NULL
        } else{
            efo.selection <- unlist(get_selected(tree))
            efo.split <- strsplit(efo.selection,':')[[1]]
            ## efo.id <- paste(efo.split[1],efo.split[2],sep='_')
            efo.id <- efo.split[1]

            ## efo.req <- dbFetch(dbSendQuery(mercator.db.con,sprintf("SELECT termtree FROM efo_tree WHERE id='%s'",efo.id)))
            ## efo.info <- fromJSON(efo.req$termtree)

            efo.info <- fromJSON(sprintf('http://localhost:3000/efo_info/%s',efo.id))
            colVar <- rep('unlabelled',length(tsne.order))
            names(colVar) <- tsne.order

            for(label in names(efo.info)){
                used.ids <- intersect(tsne.order,efo.info[[label]])
                colVar[used.ids] <- label
            }
        }
        return(colVar)
    })
    
    meshVec <- reactive({

        tree <- input$meshTree

        if(is.null(tree) | length(get_selected(tree)) == 0){
            colVar=NULL
        } else{
            mesh.selection <- unlist(get_selected(tree))
            mesh.id <- strsplit(mesh.selection,':')[[1]][1]
            mesh.info <- fromJSON(sprintf('http://localhost:3000/ontology_info/%s',mesh.id))

            ## mesh.req <- dbFetch(dbSendQuery(mercator.db.con,sprintf("SELECT termtree FROM mesh_tree WHERE id='%s'",mesh.id)))
            ## mesh.info <- fromJSON(mesh.req$termtree)
            
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

            results <- rbind(results,data.frame(samps=entry,group=rep(rowLabel,length(entry)),stringsAsFactors=F))


        }


        return(results)

    })

    projectionVec <- reactive({

        rows <- input$sampleInputTable_rows_selected

        if(is.null(rows)){
            return(NULL)
        }

        entry <- user.selections$gene.input.results[[rows]]

        return(entry)

    })

    selectionBarVec <- reactive({

        rows <- input$selectionList_rows_selected

        if(is.null(rows) & length(user.selections$selection.list) > 0){
            rows <- 1:length(user.selections$selection.list)
        }

        indVec <- c()
        groupVec <- c()

        dataResults <- data()

        results <- data.frame()

        for(rowNum in rows){

            entry <- user.selections$selection.list[[rowNum]]

            rowLabel <- sprintf('Selection: %s',user.selections$selection.datalist[rowNum,'name'])


            results <- rbind(results,data.frame(samps=entry,group=rep(rowLabel,length(entry))))


        }

        


        return(results)

    })


    violinColVar <- reactive({
        colVar <- NULL
        if(input$colorButton == 0 | input$violinXFactors %in% c('No Coloring','Gene','Mesh','Tissue','DOID','efo','KMeans','Louvain','Projection')){
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
        if(input$colorButton == 0 | input$colorfactors %in% c('No Coloring','Gene','Mesh','Tissue','DOID','efo','KMeans','Louvain','Projection')){
            return(colVar)
        }

        colVar = apply(data()[,input$colorfactors,drop=FALSE],1,paste,collapse='+')
            
        ## }

        return(colVar)
        ## }
    })

    barColVar <- reactive({
        
        colVar <- NULL
        ## if(TRUE){return(colVar)}
        if(input$colorButton == 0 | input$barPlotXaxis %in% c('No Coloring','Gene','Mesh','Tissue','DOID','efo','KMeans','Louvain','Projection','marker_tissue_type','marker_cancer_type','marker_cell_type','marker_cell_name')){
            return(colVar)
        }

        colVar = apply(data()[,input$barPlotXaxis,drop=FALSE],1,paste,collapse='+')
            
        return(colVar)

    })

    output$plotlyClick <- renderUI({

        click.event <- event_data('plotly_click',source='tsne')

        if(is.null(click.event) == T) return('Click on a t-SNE point to link to SRA entry')

        curveNumber <- click.event$curveNumber+1
        inds <- click.event$pointNumber+1

        run.id <- user.selections$tsne.traces[[curveNumber]][inds]

        samp.id <- tsne.meta[run.id,'samp_id']

        case.id <- tsne.meta[run.id,'tcga.case.ids$gdc_cases.case_id']

        if(tsne.meta[run.id,'proj_id'] == 'TCGA'){
            ## url <- a(sprintf('TCGA link to sample %s',samp.id),href=paste('https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22',samp.id,'%22%5D%7D%7D%5D%7D',sep=''),target='_blank')
            
            url <- a(sprintf('Link to TCGA case %s',samp.id),href=sprintf('https://portal.gdc.cancer.gov/cases/%s',case.id),target='_blank')
        } else{

            url <- a(sprintf('Link to SRA run %s',user.selections$tsne.traces[[curveNumber]][inds]),href=sprintf('https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=%s',user.selections$tsne.traces[[curveNumber]][inds]),target='_blank')

        }

        return(tagList(url))

        ## print(head(user.selections$tsne.traces[[curveNumber]]))

        ## print(user.selections$tsne.traces[[curveNumber]][inds])

    })

    getGeneTable <- reactive({

        if(input$geneGroup == 'all'){
            dt <- marker.list[['all']]
        } else{
            dt <- marker.list[[as.double(input$geneGroup)+1]]
        }

    })

    getSampleTable <- reactive({

        selectionRes <- selectionVec()
        louvainVec <- louvainVec()
        meshResults <- meshVec()
        tissueResults <- tissueVec()
        doidResults <- doidVec()
        efoResults <- efoVec()
        dt <- tsne.meta[,c(1,2,3,4,5,6,7)]
        colnames(dt) <- c('proj_id','samp_id','run_id','sample_type','tissue_general','tissue_detail','tcga_case_id')

        dt$Louvain <- louvainVec

        if(is.null(meshResults)){
            dt$mesh <- NA
        } else{
            dt$mesh <- meshResults
        }

        if(is.null(tissueResults)){
            dt$tissue <- NA
        } else{
            dt$tissue <- tissueResults
        }

        if(is.null(doidResults)){
            dt$doid <- NA
        } else{
            dt$doid <- doidResults
        }

        if(is.null(efoResults)){
            dt$efo <- NA
        } else{
            dt$efo <- efoResults
        }

        if(input$sampleTableControl == 'Selection'){
            dt <- dt[selectionRes$samps,]
            dt$group <- selectionRes$group
        }

        rownames(dt) <- c()

        return(dt)

    })

    output$sampleTableDownload <- downloadHandler(
        filename = function() {
            paste('samples-',Sys.Date(),'.csv',sep='')
        },
        content = function(con) {
            write.csv(getSampleTable(),con)
        }
    )

    output$geneTableDownload <- downloadHandler(
        filename = function() {
            paste('gene-table-',input$geneGroup,'-',Sys.Date(),'.csv',sep='')
        },
        content = function(con) {
            write.csv(getGeneTable(),con)
        }
    )

    mesh.tree.dat <- readRDS('data/mesh_tree_flat.RDS')
    output$meshTree <- renderTree(mesh.tree.dat)

    tissue.tree.dat <- readRDS('data/tissue_flat_tree.RDS')
    output$tissueTree <- renderTree(tissue.tree.dat)

    doid.tree.dat <- readRDS('data/doid_tree_flat.RDS')
    output$doidTree <- renderTree(doid.tree.dat)

    efo.tree.dat <- readRDS('data/efo_tree_flat_pc3sd.RDS')
    output$efoTree <- renderTree(efo.tree.dat)

    output$metadataBar <- renderPlot({
        
        input$colorButton

        isolate({
            ## euclid.file <- input$euclid_input
            ## spear.file <- input$spearman_input
            ## euclid.pca.file <- input$euclid_pca_input
            colVarResults <- barColVar()
            dataResults <- data()
            meshResults <- meshVec()
            tissueResults <- tissueVec()
            doidResults <- doidVec()
            efoResults <- efoVec()
            barCategory <- input$barPlotXaxis
            barGroup <- input$barPlotFactor
            selectionRes <- selectionBarVec()
            ## kmeansVec <- kmeansVec()
            louvainVec <- louvainVec()            
        })

        if(input$colorButton == 0 | barCategory == 'No Coloring'){
            ## plot.dat <- data.frame(Label= c('All'), number = c(length(tsne.order)),stringsAsFactors=F)            
            xGroup <- rep('All',length(tsne.order))
        } else if(barCategory == 'Louvain'){
            xGroup <- louvainVec
        } else if(barCategory == 'Mesh'){
            xGroup <- meshResults
        } else if(barCategory == 'Tissue'){
            xGroup <- tissueResults
        } else if(barCategory == 'DOID'){
            xGroup <- doidResults
        } else if(barCategory == 'efo'){
            xGroup <- efoResults
        } else if(barCategory == 'marker_tissue_type'){
            louvain.group <- as.numeric(gsub('Louvain Cluster ','',barGroup))
            xGroup <- cellmarker.annotation.counts[[louvain.group+1]][['yes']][['marker_tissue_type']]
        } else if(barCategory == 'marker_cancer_type'){
            louvain.group <- as.numeric(gsub('Louvain Cluster ','',barGroup))
            xGroup <- cellmarker.annotation.counts[[louvain.group+1]][['yes']][['marker_cancer_type']]
        } else if(barCategory == 'marker_cell_type'){
            louvain.group <- as.numeric(gsub('Louvain Cluster ','',barGroup))
            xGroup <- cellmarker.annotation.counts[[louvain.group+1]][['yes']][['marker_cell_type']]
        } else if(barCategory == 'marker_cell_name'){
            louvain.group <- as.numeric(gsub('Louvain Cluster ','',barGroup))
            xGroup <- cellmarker.annotation.counts[[louvain.group+1]][['yes']][['marker_cell_name']]
        }else{
            xGroup <- colVarResults
        }

        ## print(xGroup)

        marker.flag <- barCategory %in% c('marker_cell_name','marker_cell_type','marker_cancer_type','marker_tissue_type')

        if(input$colorButton == 0 | barGroup == 'All' | marker.flag){}
        else if(barGroup == 'Selections'){
            xGroup <- xGroup[as.character(selectionRes$samps)]
        }else if(grepl('Kmeans',barGroup)){
            kmeans.group <- as.integer(gsub('Kmeans Cluster ','',barGroup))
            used.samps <- names(which(kmeansVec==kmeans.group))
            xGroup <- xGroup[as.character(used.samps)]

        }else if(grepl('Louvain',barGroup)){
            louvain.group <- gsub('Louvain Cluster ','',barGroup)
            used.samps <- names(which(louvainVec==louvain.group))
            xGroup <- xGroup[as.character(used.samps)]
        }

        ## labels <- unlist(lapply(xGroup,function(x) strsplit(x,' [+] ')[[1]]))
        labels <- xGroup

        plot.title <- ''

        if(!marker.flag){
            metadataTable <- table(labels)
        }else{
            if(sum(labels) == 0){
                metadataTable <- c('x'=0)
                plot.title=sprintf('No unique markers for Cluster %d',louvain.group)
            } else{
                metadataTable <- labels
            }
        }

        ## metadataTable <- table(xGroup)

        plot.dat <- data.frame(Label = names(metadataTable), number = as.vector(metadataTable),stringsAsFactors=F)


        plot.dat$Label[plot.dat$Label=='NA'] <- 'unlabelled'

        ## table(plot.dat$Label)
        
        x.labels <- gsub("(.{14,}?)\\s","\\1\n",plot.dat$Label)
        names(x.labels) <- unique(plot.dat$Label)

        if('Other' %in% x.labels){
            x.labels <- x.labels[x.labels != 'Other']
            x.labels <- c(x.labels,'Other'='Other')
        } else if('unlabelled' %in% x.labels){
            x.labels <- x.labels[x.labels != 'unlabelled']
            x.labels <- c(x.labels,'unlabelled'='Unlabelled')
        }

        print(head(plot.dat))

        output.plot <- ggplot() +
            geom_bar(data=plot.dat,aes(x=Label,y=number,fill=Label),stat='identity') +
            theme(panel.background= element_blank(),
                  axis.text.x = element_text(size=14,angle=45,vjust=0.5),
                  axis.text.y = element_text(size=18),
                  legend.text = element_text(size=24),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size=20),
                  legend.title = element_blank(),
                  legend.position='none',
                  plot.title=element_text(hjust=0.9,size=22)
                  ) +
            scale_x_discrete(labels=x.labels,limits=names(x.labels)) +
            ylab('') +
            ggtitle(plot.title)

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
            tissueResults <- tissueVec()
            doidResults <- doidVec()
            efoResults <- efoVec()
            colorFactors <- input$violinXFactors
            ## kmeansVec <- kmeansVec()
            louvainVec <- louvainVec()
            selectionRes <- selectionVec()
            ## projectionVec <- user.selections$gene.input.results
            projectionVec <- projectionVec()
            yFactors <- input$violinYFactors
            gene.id <- input$whichGene
            sampleInputSelected <- input$sampleInputTable_rows_selected
            projection.list <- projection.getDatalist()
        })

        ## plot.dat <- cbind(dataResults,colVarResults)
        ## colnames(plot.dat) <- c('y1','y2','color')

        if(input$colorButton == 0 | colorFactors == 'No Coloring'){
            xGroup <- 'All'
        } else if(colorFactors == 'Mesh'){
            xGroup <- meshResults
        } else if(colorFactors == 'Tissue'){
            xGroup <- tissueResults
        } else if(colorFactors == 'DOID'){
            xGroup <- doidResults
        } else if(colorFactors == 'efo'){
            xGroup <- efoResults
        }else if(colorFactors == 'KMeans'){
            xGroup <- kmeansVec
            ## xGroup <- colVarResults
            ## xGroup <- colorFactors
        } else if(colorFactors == 'Louvain'){
            xGroup <- louvainVec
        } else{
            xGroup <- colVarResults
        }
        ## if(is.null(xGroup)){
        ##     yGroup <- NULL
        ## }
        ## else{
        
        if(yFactors == 'projection'){
            yGroup <- -1 * projectionVec
            ylabel='Projection'
            ## print(sampleInputSelected)
            ## print(projection.list[
            plot.title=paste('Projection of ',projection.list[sampleInputSelected,'name'],sep='')
        } else{
            ## yGroup <- log10(geneVecResults+1)
            yGroup <- geneVecResults
            gene.label <- names(gene.choices)[gene.num.vec[gene.id]]
            ylabel='log2( TPM + 1 )'
            ## ylabel <- paste('log10( TPM + 1) of ',gene.label,sep='')
            plot.title=paste('Expression of ',gene.label,sep='')
        }

        
        if(any(xGroup == 'All')){
            ## used.names <- names(yGroup)
            plot.dat <- data.frame(x=as.factor(xGroup),y=yGroup)
            sizeVec <- rep(3,nrow(plot.dat))
            top.groups <- NULL
        } else if(colorFactors == 'Louvain'){

            ## if(xGroup=='All'){
            ##     used.names <- names(yGroup)
            ## } else{
            used.names <- intersect(names(xGroup),names(yGroup))
            ## }

            ## print(top.clus.per.genes[[gene.id]])

            plot.dat <- data.frame(x=as.character(xGroup[used.names]),y=yGroup[used.names],stringsAsFactors=F)

            group.means <- aggregate(plot.dat[,'y'],list(plot.dat[,'x']),mean)
            group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]
            ## print(group.means[1:10])

            ## print(top.clus.per.genes[[gene.id]])
            
            if(!is.null(selectionRes)){

                selectionRes <- selectionRes[selectionRes$samps %in% used.names,]

                num.selections <- length(unique(selectionRes$group))

                ## plot.dat <- rbind(plot.dat,data.frame(x=as.factor(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))

            } else{
                selectionRes <- data.frame()
                num.selections <- 0
            }


            top.groups <- group.means[1:(10-num.selections)]
            ## top.groups <- top.clus.per.genes[[gene.id]][1:(length(top.clus.per.genes[[gene.id]])-num.selections)]-1

            ## plot.dat[which(! (plot.dat$x %in% as.character(top.clus.per.genes[[gene.id]][1:(length(top.clus.per.genes[[gene.id]])-num.selections)-1,'x'] <- 'Other'
            plot.dat[which(! (plot.dat$x %in% as.character(top.groups))),'x'] <- 'Other'

            plot.dat <- rbind(plot.dat,data.frame(x=as.factor(selectionRes$group),y=yGroup[selectionRes$samps]))

            sizeVec <- rep(3,nrow(plot.dat))
            sizeVec[which(plot.dat$x != 'Other')] <- 5

            top.groups <- c(top.groups,unique(selectionRes$group),'Other')

        } else if(colorFactors=='tissue_general'){

            used.names <- intersect(names(xGroup),names(yGroup))

            plot.dat <- data.frame(x=as.character(xGroup[used.names]),y=yGroup[used.names],stringsAsFactors=F)

            plot.dat <- subset(plot.dat,x != 'NA')

            if(!is.null(selectionRes)){
                selectionRes <- selectionRes[selectionRes$samps %in% used.names,]
            }

            plot.dat <- rbind(plot.dat,data.frame(x=as.factor(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))

            sizeVec <- rep(3,nrow(plot.dat))

            top.groups <- NULL

        } else if(colorFactors == 'tissue_detail'){
            
            used.names <- intersect(names(xGroup),names(yGroup))
            plot.dat <- data.frame(x=as.character(xGroup[used.names]),y=yGroup[used.names],stringsAsFactors=F)

            plot.dat <- subset(plot.dat,x != 'NA')

            if(!is.null(selectionRes)){

                selectionRes <- selectionRes[selectionRes$samps %in% used.names,]

                num.selections <- length(unique(selectionRes$group))

                ## plot.dat <- rbind(plot.dat,data.frame(x=as.factor(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))

            } else{

                selectionRes <- data.frame()
                num.selections <- 0
            }

            sizeVec <- rep(3,nrow(plot.dat))

            group.means <- aggregate(plot.dat[,'y'],list(plot.dat[,'x']),mean)
            group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]
            group.means <- group.means[! group.means %in% c('Skin','Soft Tissue','Thymus')]

            top.groups <- group.means[1:(13 - num.selections)]

            plot.dat[which(! (plot.dat$x %in% as.character(top.groups))),'x'] <- 'Other'

            plot.dat <- rbind(plot.dat,data.frame(x=as.character(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))

            top.groups <- c(top.groups,unique(selectionRes$group),'Other')
            
            sizeVec[which(plot.dat$x != 'Other')] <- 5
            

        }
        ## else if(any(xGroup != 'All')){
        else{
            used.names <- intersect(names(xGroup),names(yGroup))
            plot.dat <- data.frame(x=as.character(xGroup[used.names]),y=yGroup[used.names],stringsAsFactors=F)

            labelled.dat <- subset(plot.dat,x != 'unlabelled')

            group.cnts <- table(labelled.dat[,'x'])
            used.groups <- names(which(group.cnts>(nrow(labelled.dat) / 100)))
            labelled.dat <- subset(labelled.dat,x %in% used.groups)


            if(!is.null(selectionRes)){

                selectionRes <- selectionRes[selectionRes$samps %in% used.names,]

                num.selections <- length(unique(selectionRes$group))

                ## plot.dat <- rbind(plot.dat,data.frame(x=as.factor(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))

            } else{

                selectionRes <- data.frame()
                num.selections <- 0
            }

            sizeVec <- rep(3,nrow(plot.dat))

            top.groups <- NULL

            if((length(unique(xGroup)) + num.selections) > 10){

                group.means <- aggregate(labelled.dat[,'y'],list(labelled.dat[,'x']),mean)
                group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]

                ## group.means <- group.means[group.means %in% used.groups]

                top.groups <- group.means[1:(10 - num.selections)]

                plot.dat[which(! (plot.dat$x %in% as.character(top.groups))),'x'] <- 'Other'

                sizeVec[which(plot.dat$x != 'Other')] <- 5
                
                top.groups <- c(top.groups,unique(selectionRes$group),'Other')

            }

            

            plot.dat <- rbind(plot.dat,data.frame(x=as.character(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))
            
            sizeVec[which(plot.dat$x != 'Other' & plot.dat$x != 'unlabelled')] <- 5

        } 




        x.labels <- gsub("(.{14,}?)\\s","\\1\n",unique(plot.dat$x))
        names(x.labels) <- unique(plot.dat$x)

        if('Other' %in% x.labels){
            x.labels <- x.labels[x.labels != 'Other']
            x.labels <- c(x.labels,'Other'='Other')
        } else if('unlabelled' %in% x.labels){
            x.labels <- x.labels[x.labels != 'unlabelled']
            x.labels <- c(x.labels,'unlabelled'='Unlabelled')
        } else if('Eye' %in% x.labels){
            x.labels <- x.labels[x.labels != 'Eye']
            x.labels <- c(x.labels,'Eye'='Eye')
        }

        if(!is.null(top.groups)){

            x.labels <- x.labels[top.groups]

        }

        ## min.y <- min(subset(plot.dat,y>0)$y)

        output.plot <- ggplot() +
            geom_violin(data=plot.dat,aes(x=x,y=log2(abs(y)+1),fill=x),scale='width') +
            geom_jitter(data=plot.dat,aes(x=x,y=log2(abs(y)+1)),width=0.1,alpha=0.1,colour='gray',size=sizeVec) +
            theme(##panel.background= element_blank(),
                axis.text.x = element_text(size=14,angle=45,vjust=0.5),
                axis.text.y = element_text(size=18),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size=20),
                legend.title = element_blank(),
                legend.position='none',
                plot.title=element_text(hjust=1.0,size=22)
                  ) +
            ylab(ylabel) +
            ggtitle(plot.title) +
            scale_x_discrete(labels=x.labels,limits=names(x.labels))

            ## scale_
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
            tissueResults <- tissueVec()
            doidResults <- doidVec()
            efoResults <- efoVec()
            colorFactors <- input$colorfactors
            ## kmeansVec <- kmeansVec()
            louvainVec <- louvainVec()
            ## projectionVec <- user.selections$gene.input.results
            projectionVec <- projectionVec()
            ## input$tree
            ## input$colorfactors
            markerSize <- input$markerSize
        })


        ## print(colorFactors)

        ## color.ramp <- c('#800080','#FFFF00')

        color.ramp <- NULL
        
        if(input$colorButton == 0 | colorFactors == 'No Coloring'){
            colVarPlot <- NULL

        } else if(colorFactors == 'Mesh'){
            colVarPlot <- meshResults

        } else if(colorFactors == 'Tissue'){
            colVarPlot <- tissueResults
            
        } else if(colorFactors == 'DOID'){
            colVarPlot <- doidResults

        } else if(colorFactors == 'efo'){
            colVarPlot <- efoResults

        } else if(colorFactors == 'Gene'){
            colVarPlot <- log2(geneVecResults+1)

        } else if(colorFactors == 'KMeans'){
            colVarPlot <- as.factor(kmeansVec[rownames(dataResults)])

        } else if(colorFactors == 'Louvain'){
            colVarPlot <- as.factor(louvainVec[rownames(dataResults)])

        } else if(colorFactors == 'Projection'){
            colVarPlot <- log2(projectionVec+1)

            ## bottom.ramp <- colorRampPalette(c('#D73027','#FFFFBF'))
            ## bottom.ramp <- colorRampPalette(c('#FF0000','#FFFFBF'))
            bottom.ramp <- colorRampPalette(c('#FF0000','#ADD8E6'))

            ## top.ramp <- colorRampPalette(c('#FFFFBF','#4575B4'))
            ## top.ramp <- colorRampPalette(c('#FFFFBF','#FFFFF8'))
            top.ramp <- colorRampPalette(c('#ADD8E6','#F6FBFC'))

            ## min.point <- min(min(colVarPlot),0.1439294) ### using true minimum, is waaaaaay too low
            min.point <- min(min(colVarPlot),13.51984) # 1st percentile
            ## min.point <- min(min(colVarPlot),10.85842) # testing, may be too small
            ## max.point <- max(max(colVarPlot),25.65136)
            max.point <- 25.65136
            
            ## mean.point <- min(mean(colVarPlot),16.69396)
            ## mean.point <- 16.69396
            mean.point <- mean(colVarPlot)

            ## print(min.point)
            ## print(max.point)

            mid.point <- round((mean.point - min.point) / (max.point - min.point)*100)

            ## print(mid.point)

            ## mid.point <- round((mean(colVarPlot) - min(colVarPlot)) / (max(colVarPlot) - min(colVarPlot)) * 100)

            color.ramp <- c(bottom.ramp(mid.point),top.ramp(100-mid.point))

            sample.min.point <- round((min(colVarPlot) - min.point) / (max.point - min.point)*100)

            ## print(sample.min.point)

            color.ramp <- color.ramp[sample.min.point:100]
            
        }
        else{
            colVarPlot <- colVarResults
        }

        ## dataResults$run.id <- rownames(dataResults)
        ## dataResults <- cbind(dataResults,KMeans = as.character(kmeansVec[rownames(dataResults)]))
        dataResults <- cbind(dataResults,Louvain = as.character(louvainVec[rownames(dataResults)]))

        if(length(geneVecResults) > 0){
            dataResults <- cbind(dataResults,Gene = as.character(geneVecResults))
        } else{
            dataResults <- cbind(dataResults,Gene=NA)
        }

        colVarAnnotation <- colVarPlot

        if(typeof(colVarPlot) == 'character'){
            colVarPlot <- gsub("(.{14,}?)\\s","\\1<br>",colVarPlot)
        }

        ## if(!is.null(colVarPlot)){

        if(any(colVarPlot == 'NA') | any(colVarPlot == 'unlabelled')){

            colVarPlot[colVarPlot == 'NA'] <- 'unlabelled'
            colVarAnnotation[colVarAnnotation == 'NA'] <- 'unlabelled'

            set.seed(11)
            
            color.ramp <- setNames(colorRampPalette(c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F"))(length(unique(colVarPlot))),
                                   ## colorRampPalette(brewer.pal(12,"Set3"))(length(unique(colVarPlot))),
                                   sample(unique(colVarPlot)))
            color.ramp['unlabelled'] <- '#F0F0F0'
            
        }

        if(!is.null(colVarPlot) & colorFactors != 'Gene' & colorFactors != 'Projection'){
            
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

        output.plot <- plot_ly(dataResults, x = ~y1, y = ~y2,mode="markers",type='scattergl',color = colVarPlot,
                               colors=color.ramp,
                               hoverinfo='text',
                               ## opacity=plot.opacity,
                               ## group = colVarPlot,
                               text=~paste('ID: ',run_id,
                                           '<br>Project: ',proj_id,
                                           '<br>Tissue General: ',tissue_general,
                                           '<br>Tissue Detail: ',tissue_detail,
                                           '<br>Sample Type: ',sample_type,
                                           '<br>Sample ID: ',samp_id,
                                           ## '<br>Kmeans: ',KMeans,
                                           '<br>Louvain: ',Louvain,
                                           '<br>Gene: ',Gene,
                                           '<br>Annotation: ',colVarAnnotation),
                               marker = list(size=markerSize),
                               source='tsne') %>%
            ## config(p = .,modeBarButtonsToRemove = c("zoom2d",'toImage','autoScale2d','hoverClosestGl2d'),collaborate=FALSE,cloud=FALSE) %>%
            ## config(collaborate=FALSE) %>%
            layout(dragmode = "pan",
                   xaxis=ax,
                   yaxis=ax,
                   legend=list(font=list(family='sans-serif',size=11))) ## %>%

        output.plot
        

    })
  
})

