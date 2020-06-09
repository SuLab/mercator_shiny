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

    ## tsne.data <- readRDS('../../results/recount/tsne/recount_tsne_list_pca_noNAs_pcinit_eta200_over50_pc3sd_poscounts_log.RDS')[['30.50']]

    ## tsne.data <- readRDS('data/tsne_recount_pc3sd_poscounts_p175_nv750.RDS')
    ## tsne.data <- readRDS('../../results/recount/tmp/tsne_over50_poscounts_log_combp_30_250.RDS')

    ## tsne.data <- readRDS('data/recount_tsne_pca_noNAs_pcinit_eta2875_over50_pc3sd_tpm_log_p60_nv750.RDS')

    ## tsne.data <- readRDS('../../results/recount/tsne/recount_tsne_list_pca_pcinit_eta2875_over50_pc3sd_tpm_log_90th_var_genes.RDS')[['40.100']]
    ## tsne.data <- readRDS('../../results/recount/tsne/recount_tsne_list_pca_noNAs_pcinit_eta2875_over50_pc3sd_tpm_log.RDS')[['60.100']]
    tsne.data <- readRDS('data/recount_tsne_pca_pcinit_eta2875_over50_pc3sd_tpm_log_90th_var_genes_p40_nv100.RDS')
    
    ## tsne.data <- tsne.data + matrix(data=rnorm(2*nrow(tsne.data),sd=0.5),ncol=2)

    tsne.order <- rownames(tsne.data)

    ## tsne.meta <- readRDS('data/recount2_meta_edited_hist_pc3sd_over50.RDS')
    tsne.meta <- readRDS('data/recount2_meta_edited_hover_labels.RDS')

    ## tsne.meta <- readRDS('data/recount_meta_sampletype_urls.RDS')
    ## rownames(tsne.meta) <- tsne.meta$run_id
    ## tsne.meta <- tsne.meta[tsne.order,]

    ## tcga.gtex.meta <- readRDS('data/tcga_gtex_meta_edited_histologies.RDS')

    ## rownames(tcga.gtex.meta) <- tcga.gtex.meta$data_id
    ## tcga.gtex.meta <- tcga.gtex.meta[tsne.order,]

    ## tsne.meta <- cbind(tsne.meta,tcga.gtex.meta[,c('tissue_general','tissue_detail')])

    gene.tpm.samps <- readRDS('data/mrmnorm_mat_rownames.RDS')

    gene.choices <- readRDS('data/shiny_gene_choice_list.RDS')
    gene.num.vec <- readRDS('data/shiny_gene_num_vec.RDS')

    ## louvain.dat <- read.table('../../results/recount/clustering/leiden_rsweep_pca_over50_pc3sd_tpm_log_nv100_narrow_90th_var_genes_k30_sim.tsv',sep='\t',row.names=1)
    ## ## colnames(louvain.dat) <- c('0.1', '0.075', '0.05', '0.025', '0.01', '0.0075', '0.005', '0.0025', '0.001', '0.00075', '0.0005', '0.00025', '0.00001')
    ## ## colnames(louvain.dat) <- c('0.1','0.11','0.15','0.175','0.2')
    ## colnames(louvain.dat) <- c('0.007', '0.0075', '0.008', '0.0085', '0.009', '0.0095', '0.01', '0.0105', '0.011', '0.0115', '0.012', '0.0125', '0.013')
    ## louvain.dat <- read.table('../../results/recount/clustering/leiden_rsweep_narrow_pca_over50_pc3sd_tpm_log_k20_sim.tsv',sep='\t',row.names=1)
    ## louvain.dat <- read.table('../../results/recount/clustering/leiden_rsweep_narrow_pca_over50_pc3sd_tpm_log_nv100_k30_sim.tsv',sep='\t',row.names=1)
    ## ## colnames(louvain.dat) <- c('0.005','0.0045','0.004','0.0035','0.003','0.0025')
    ## colnames(louvain.dat) <- c('0.03','0.0275','0.025','0.0225','0.02','0.0175','0.015','0.0125','0.01')
    ## colnames(louvain.dat) <- c('0.1', '0.075', '0.05', '0.025', '0.01', '0.0075', '0.005', '0.0025', '0.001', '0.00075', '0.0005', '0.00025', '0.00001') 

    ## louvain.vec <- louvain.dat[,'0.009']
    ## names(louvain.vec) <- tsne.order
    louvain.vec <- readRDS('data/leiden_r9e-3_over50_pc3sd_tpm_log_k30_sim_90th_var_genes.RDS')
    ## louvain.vec <- readRDS('data/leiden_r25e-3_over50_pc3sd_poscounts_k40_sim_nosingles.RDS')
    ## louvain.vec <- readRDS('data/leiden_pca_r5e-3_pc3sd_tpm_log_k20_sim_nosingles.RDS')
    ## louvain.vec <- readRDS('data/leiden_r25e-3_pc3sd_tpm_log_k40_sim_nosingles.RDS')
    ## louvain.vec <- readRDS('data/leiden_r9e-3_over50_pc3sd_tpm_log_k30_sim_90th_var_genes.RDS')
    louvain.choices <- sort(unique(louvain.vec))
    louvain.choices <- sapply(louvain.choices,function(x) sprintf('Louvain Cluster %s',x))
    names(louvain.choices) <- louvain.choices

    cellmarker.info <- readRDS('data/cell_marker_gene_table_fixed_ens.RDS')
    rownames(cellmarker.info) <- cellmarker.info$V3

    cellmarker.annotation.counts <- readRDS('data/cell_marker_annotations_marker_counts.RDS')

    ## wgcna.dat <- readRDS('data/wgcna_app_table_pos_log_pc3sd_min100_ds2.RDS')

    ## onevent('item_remove','hooHa-selectized',print('bar'))

    user.selections <- reactiveValues(
        selection.list = list(),
        selection.datalist = data.frame(),
        gene.input.results = list(),
        gene.input.datalist = data.frame(),
        tsne.traces = list(tsne.order),
        previous.tab = 't-SNE',
        violinX.change = TRUE,
        colorButtonActive = FALSE,
        currentViolinGroups = c(),
        violinGroupRendered = FALSE,
        firstViolinRender = FALSE
    )

    updateSelectizeInput(session,
                         'whichGene',
                         choices=gene.choices,
                         server=TRUE,
                         options=list(maxOptions=10,
                                      closeAfterSelect=TRUE
                                      )
                         )


    onevent('change','violinXFactors',function(event){
        
        isolate({user.selections$violinX.change <- TRUE})

    })
    ## updateSelectizeInput(session,
    ##                      'hooHa',
    ##                      choices=c(1,2,3,4,5),
    ##                      selected=c(1,2),
    ##                      server=TRUE,
    ##                      options=list(maxOptions=10,
    ##                                   closeAfterSelect=FALSE)
    ##                      )


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
                              label='Bar Plot X-axis',                              
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
                              label='Bar Plot X-axis',
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


    observeEvent(input$wgcnaGeneTable_rows_selected, {

        selected.gene <- getWgcnaGeneTable()[input$wgcnaGeneTable_rows_selected,3]

        updateSelectizeInput(session,
                             'whichGene',
                             choices = gene.choices,
                             server=T,
                             options=list(maxOptions=10,
                                          closeAfterSelect=T
                                          ),
                             selected=gene.choices[gene.num.vec[selected.gene]]
                             )

    })

    observeEvent(input$geneTable_rows_selected, {

        ## if(input$geneGroup == 'all'){
        ##     ## selected.marker <- rownames(marker.list[['all']])[input$geneTable_rows_selected]
        ## }
        ## else{
        ##     selected.marker <- rownames(marker.list[[as.double(input$geneGroup)+1]])[input$geneTable_rows_selected]
        ## }

        selected.marker <- rownames(getGeneTable())[input$geneTable_rows_selected]

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

        ## saveRDS(event.data,'~/Documents/Projects/Mercator/tmp/event_dat.RDS')

        if(is.null(event.data) == TRUE) return(NULL)

        ## curveNumber <- event.data$curveNumber+1

        ## inds <- event.data$pointNumber+1

        selection.vec <- c()

        for(curveNum in unique(event.data$curveNumber)){

            inds <- subset(event.data,curveNumber==curveNum)$pointNumber + 1
            
            selection.vec <- c(selection.vec,user.selections$tsne.traces[[curveNum+1]][inds])

        }

        ## user.selections$selection.list[[length(user.selections$selection.list)+1]] <- user.selections$tsne.traces[[curveNumber[1]]][inds]
        user.selections$selection.list[[length(user.selections$selection.list)+1]] <- selection.vec

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

    output$wgcnaGeneTable <- DT::renderDT({

        wgcna.tab <- getWgcnaGeneTable()

        ## dataDT <- DT::datatable(wgcna.tab[,c(11,9,10,4,5,6,7,8)],
        dataDT <- DT::datatable(wgcna.tab[,c(11,9,10,4)],
                                selection='single',
                                options=list(
                                    'searching'=TRUE,
                                    'legnthChange'=FALSE,
                                    'pagingType'='full',
                                    'autoWidth'=TRUE,
                                    'pageLength'=25),
                                escape=FALSE,
                                rownames=FALSE
                                )

        return(dataDT)
    })
                                                                
    output$wgcnaGoTable <- DT::renderDT({

        go.tab <- getWgcnaGoTable()

        dataDT <- DT::datatable(go.tab,
                                selection='none',
                                options=list(
                                    'searching'=TRUE,
                                    'legnthChange'=FALSE,
                                    'pagingType'='full',
                                    'autoWidth'=TRUE,
                                    'pageLength'=25),
                                escape=FALSE,
                                rownames=FALSE
                                )

        return(dataDT)

    })

    output$geneTable <- DT::renderDT({

        ## if(input$geneGroup == 'all'){
        ##     label.1 <- 'all'
        ##     label.2 <- input$geneGroupSecond
        ## }
        ## else{
        ##     label.1 <- as.numeric(input$geneGroup)
        ##     label.2 <- as.numeric(input$geneGroupSecond)

        ##     if(label.1 > label.2){
        ##         tmp <- label.1
        ##         label.1 <- label.2
        ##         label.2 <- tmp
        ##     }
        ## }

        ## ## label.2 <- as.numeric(input$geneGroupSecond)+1

        ## label <- sprintf('%s.%s',label.1,label.2)

        ## marker.tab <- fromJSON(sprintf('http://localhost:3000/pairwise_markers/%s',label))

        ## rownames(marker.tab) <- marker.tab[,1]
        ## colnames(marker.tab) <- c('ens_id','p-val','unique-fcs','total-fcs')

        ## marker.tab <- cbind(cellmarker.info[rownames(marker.tab),],marker.tab[,c(2,3,4)])

        ## marker.tab <- marker.tab[,c(-3)]
        ## colnames(marker.tab) <- c('Symbol','ID','gene-type','tissueType','cancerType','cellType','cellName','p-val','unique-fcs','total-fcs')

        marker.tab <- getGeneTable()

        if(nrow(marker.tab) > 0){
            ## marker.tab <- marker.tab[,c(10,8,9,3,4,5,6,7,11:ncol(marker.tab))]
            marker.tab <- marker.tab[,c(10,8,9,3,11:ncol(marker.tab))]

        }

        ## print(head(marker.tab))

        ## if(input$geneGroup == 'all'){
        ##     dt <- marker.list[['all']]
        ## } else{
        ##     dt <- marker.list[[as.double(input$geneGroup)+1]]
        ## }

        ## dataDT <- DT::datatable(dt,
        dataDT <- DT::datatable(marker.tab,
                                ## colnames=c('Ensembl ID','Symbol','Entrez ID','Gene Type','CM tissueType','CM cancerType','CM cellType','CM cellName','Unique','cluster p-val','cluster fc','cluster-all fc'),
                                selection='single',
                                options=list(
                                    'searching'=TRUE,
                                    'legnthChange'=FALSE,
                                    'pagingType'='full',
                                    'autoWidth'=TRUE,
                                    'pageLength'=25),
                                escape=FALSE,
                                rownames=FALSE
                                )
        

        return(dataDT)
    })

    output$sampleTable <- DT::renderDT({

        ## input$colorButton

        ## isolate({

        ## selectionRes <- selectionVec()
        ## louvainVec <- louvainVec()
        ## meshResults <- meshVec()
        ## tissueResults <- tissueVec()
        ## doidResults <- doidVec()
        ## efoResults <- efoVec()

        ## ## })

        ## dt <- tsne.meta[,c(1,2,3,4,5,7,8)]

        ## dt$Louvain <- louvainVec

        ## if(is.null(meshResults)){
        ##     dt$mesh <- NA
        ## } else{
        ##     dt$mesh <- meshResults
        ## }

        ## if(is.null(tissueResults)){
        ##     dt$tissue <- NA
        ## } else{
        ##     dt$tissue <- tissueResults
        ## }

        ## if(is.null(doidResults)){
        ##     dt$doid <- NA
        ## } else{
        ##     dt$doid <- doidResults
        ## }

        ## if(is.null(efoResults)){
        ##     dt$efo <- NA
        ## } else{
        ##     dt$efo <- efoResults
        ## }

        ## if(input$sampleTableControl == 'Selection'){
        ##     dt <- dt[selectionRes$samps,]
        ##     dt$group <- selectionRes$group
        ## ## } else if('Louvain' %in% input$sampleTableControl){
        ## } else if(grepl('Louvain',input$sampleTableControl)){
        ##     louvain.id <- gsub('Louvain Cluster ','',input$sampleTableControl)

        ##     dt <- subset(dt,Louvain==as.numeric(louvain.id))
        ## }

        ## rownames(dt) <- c()



        ## dataDT <- DT::datatable(dt,
        dataDT <- DT::datatable(getSampleTable()[,c(-6)],
                                options = list(pageLength=25
                                               ),
                                escape=FALSE
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

        if(input$topNavBar == 'violinPanel'){user.selections$violinGroupRendered <- TRUE}

        if(input$topNavBar == 'genePanel' | input$topNavBar == 'samplePanel' | input$topNavBar == 'wgcnaPanel'){
            
            runjs(" var $this = $('#controlPanelClickable');
                    if(!$this.hasClass('panel-collapsed')) {
                        $this.parents('.panel').find('.panel-body').slideUp();
		        $this.addClass('panel-collapsed');
		        $this.find('i').removeClass('glyphicon-chevron-up').addClass('glyphicon-chevron-down');
                    }")

        }
        else if(user.selections$previous.tab == 'genePanel' | user.selections$previous.tab == 'samplePanel' | user.selections$previous.tab == 'wgcnaPanel'){

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

            ## eu.dist <- sum((x - gene.vec)^2)

            eu.dist <- cor(x,gene.vec[1,],method='pearson')
            ## eu.dist <- cor(x,gene.vec,method='spearman')

            ## eu.dist <- (x %*% gene.vec) / (norm(x,'2') * norm(gene.vec,'2'))

            return(eu.dist)

        }

        progress <- shiny::Progress$new()

        on.exit(progress$close())

        progress$set(message='Calculating...', value=0.33)

        ## map.dat <- readRDS('data/recount_750_dim_noProj_over50_pc3sd_poscounts.RDS')
        ## rotation.dat <- readRDS('data/recount_all_big_irlba_750_over50_pc3sd_poscounts_rotmat.RDS')
        ## map.centers <- readRDS('data/recount_over50_bulkOnly_pc3sd_poscounts_colparams.RDS')
        ## gene.means <- readRDS('data/recount_over50_pc3sd_geoMeansNZ.RDS')

        ## map.dat <- readRDS('data/recount_750_dim_noProj_over50_pc3sd_tpm_log.RDS')
        ## rotation.dat <- readRDS('data/recount_all_big_irlba_750_over50_pc3sd_tpm_log_rotmat.RDS')
        ## map.centers <- readRDS('data/recount_over50_pc3sd_colparams_tpm_log.RDS')
        map.dat <- readRDS('data/recount_100_dim_noProj_over50_pc3sd_tpm_log_90th_var_genes.RDS')
        rotation.dat <- readRDS('data/recount_all_big_irlba_100_over50_pc3sd_tpm_log_90th_var_genes_rotmat.RDS')
        map.centers <- readRDS('data/recount_over50_pc3sd_colparams_tpm_log_90th_var_genes.RDS')
        gene.lens <- read.table('data/ensembl_gene_lengths.tsv',sep='\t',stringsAsFactors=F,row.names=1)

        map.centers <- subset(map.centers,vars>0)
        gene.lens <- gene.lens[rownames(map.centers),]

        gene.dat <- read.table(input$gene.vec$datapath,sep='\t',header=T)
        gene.dat <- gene.dat[rownames(map.centers),,drop=F]
        input.names <- colnames(gene.dat)

        gene.dat <- gene.dat / gene.lens
        normFactors <- colSums(gene.dat)
        gene.dat <- 1e6 * data.matrix(gene.dat) %*% diag(1/normFactors,nrow=length(normFactors))
        gene.dat <- log(gene.dat+1)
        colnames(gene.dat) <- input.names

        dist.mat <- matrix(0,nrow=nrow(map.dat),ncol=ncol(gene.dat))
        rownames(dist.mat) <- rownames(map.dat)
        colnames(dist.mat) <- colnames(gene.dat)

        for(col in colnames(gene.dat)){

            ## gene.vec <- gene.dat[[col]]
            gene.vec <- gene.dat[,col,drop=F]
            names(gene.vec) <- rownames(gene.dat)

            ## gene.ratios <- gene.vec / gene.means
            ## size.factor <- median(gene.ratios[is.finite(gene.ratios) & gene.ratios > 0])
            ## gene.vec <- gene.vec / size.factor

            center.vec <- gene.vec - map.centers$means
            center.vec <- center.vec / map.centers$vars

            rot.vec <- t(center.vec) %*% rotation.dat

            dist.vec <- apply(map.dat,1,getDistance,gene.vec=rot.vec)

            dist.mat[,col] <- dist.vec

        }

        progress$set(message='Calculated', value=0.67)

        dat.rows <- rownames(data())

        dist.mat <- dist.mat[dat.rows,,drop=F]

        for(col in colnames(dist.mat)){

            dist.vec <- dist.mat[,col]
            names(dist.vec) <- rownames(dist.mat)

            user.selections$gene.input.results[[length(user.selections$gene.input.results)+1]] <- dist.vec
            ## user.selections$gene.input.datalist <- rbind(user.selections$gene.input.datalist,data.frame('name'=col,'min'=min(log2(dist.vec+1))))
            user.selections$gene.input.datalist <- rbind(user.selections$gene.input.datalist,data.frame('name'=col,'max'=max(dist.vec)))
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

        if(input$whichGene == ''){
        ## if(TRUE){
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

        ## if(input$geneScale == 'log2gene'){
        ##     geneVec <- log2(geneVec)
        ## } else if(input$geneScale == 'log2gene1'){
        ##     geneVec <- log2(geneVec + 1)
        ## } 

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
            efo.id <- paste(efo.split[1],efo.split[2],sep='_')

            ## efo.id <- efo.split[1]

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

        ## rows <- input$selectionList_rows_selected

        ## if(is.null(rows)){
        ##     return(NULL)
        ## }

        if(length(user.selections$selection.list) == 0){

            return(NULL)

        }

        ## indVec <- c()
        ## groupVec <- c()

        ## dataResults <- data()

        results <- data.frame()

        ## for(rowNum in rows){
        for(rowNum in 1:length(user.selections$selection.list)){

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

        ## case.id <- tsne.meta[run.id,'tcga.case.ids$gdc_cases.case_id']
        case.id <- tsne.meta[run.id,'gdc_cases.case_id']

        if(is.na(run.id)) return('Click on a t-SNE point to link to SRA entry')

        if(tsne.meta[run.id,'proj_id'] == 'TCGA'){
            ## url <- a(sprintf('TCGA link to sample %s',samp.id),href=paste('https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22',samp.id,'%22%5D%7D%7D%5D%7D',sep=''),target='_blank')

            ## a(sprintf('Link to TCGA case %s',samp.id),href=sprintf('https://portal.gdc.cancer.gov/cases/%s',case.id),target='_blank')

            url <- a(sprintf('Link to TCGA case %s',samp.id),href=sprintf('https://portal.gdc.cancer.gov/cases/%s',case.id),target='_blank')

            ## print(head(tsne.meta))
            ## url <- a(sprintf('Link to TCGA case %s',user.selections$tsne.traces[[curveNumber]][inds]),href=sprintf('https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=%s',user.selections$tsne.traces[[curveNumber]][inds]),target='_blank')
        } else{

            url <- a(sprintf('Link to SRA run %s',user.selections$tsne.traces[[curveNumber]][inds]),href=sprintf('https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=%s',user.selections$tsne.traces[[curveNumber]][inds]),target='_blank')

        }

        return(tagList(url))

        ## print(head(user.selections$tsne.traces[[curveNumber]]))

        ## print(user.selections$tsne.traces[[curveNumber]][inds])

    })

    getWgcnaGoTable <- reactive({

        go.tab <- wgcna.dat[['cats']][[input$wgcnaTabControl]]

        colnames(go.tab) <- c('Category','Term','Ontology','p-val')

        return(go.tab)

    })

    getWgcnaGeneTable <- reactive({

        gene.tab <- wgcna.dat[['genes']][[input$wgcnaTabControl]]

        colnames(gene.tab) <- c('Symbol','ID','Ensembl','Gene Type','CM tissueType','CM cancerType','CM cellType','CM cellName','Gene Symbol','Entrez ID','Ensembl ID')

        return(gene.tab)
        ## return(gene.tab[,c(11,9,10,4,5,6,7,8)])

    })

    getGeneTable <- reactive({

        if(input$geneGroupSecond == 'all'){
            label.1 <- 'all'
            label.2 <- input$geneGroup

            label <- sprintf('%s.%s',label.1,label.2)

            marker.tab <- fromJSON(sprintf('http://localhost:3000/pairwise_markers/%s',label))

            rownames(marker.tab) <- marker.tab[,1]
            colnames(marker.tab) <- c('ens_id','unique','p-val','unique-fcs','total-fcs','stat')

            marker.tab <- cbind(cellmarker.info[rownames(marker.tab),],marker.tab[,c(2,3,4,5,6)],stringsAsFactors=F)

            marker.tab <- marker.tab[,c(-3)]
            colnames(marker.tab) <- c('Symbol','ID','Gene Type','CM tissueType','CM cancerType','CM cellType','CM cellName','Gene Symbol','Entrez ID','Ensembl ID','unique','p-val','unique-fcs','total-fcs','stat')

            marker.tab[['unique-fcs']] <- as.numeric(marker.tab[['unique-fcs']])
            marker.tab[['total-fcs']] <- as.numeric(marker.tab[['total-fcs']])

            ## marker.tab <- marker.tab[order(marker.tab[['unique-fcs']],decreasing=T),]
        }
        else if(input$geneGroup == input$geneGroupSecond){
            return(subset(data.frame('name'=c(1,2),'min'=c(2,3)),name>3)) ## TODO modify to have all columns
        }
        else{
            label.1 <- as.numeric(input$geneGroup)
            label.2 <- as.numeric(input$geneGroupSecond)

            mult <- -1

            if(label.1 > label.2){
                tmp <- label.1
                label.1 <- label.2
                label.2 <- tmp

                mult <- 1
            }

            label <- sprintf('%s.%s',label.1,label.2)

            marker.tab <- fromJSON(sprintf('http://localhost:3000/pairwise_markers/%s',label))

            rownames(marker.tab) <- marker.tab[,1]
            colnames(marker.tab) <- c('ens_id','p-val','fc','test.stat')
            ## colnames(marker.tab) <- c('ens_id','p-val','fc','test.stat','wilcox')

            ## marker.tab <- cbind(cellmarker.info[rownames(marker.tab),],marker.tab[,c(2,3,4,5)],stringsAsFactors=F)
            marker.tab <- cbind(cellmarker.info[rownames(marker.tab),],marker.tab[,c(2,3,4)],stringsAsFactors=F)

            marker.tab <- marker.tab[,c(-3)]
            ## colnames(marker.tab) <- c('Symbol','ID','Gene Type','CM tissueType','CM cancerType','CM cellType','CM cellName','Gene Symbol','Entrez ID','Ensembl ID','p-val','log2fc','test.stat','scaled wilcox stat')
            colnames(marker.tab) <- c('Symbol','ID','Gene Type','CM tissueType','CM cancerType','CM cellType','CM cellName','Gene Symbol','Entrez ID','Ensembl ID','p-val','log2fc','test.stat')

            marker.tab[['log2fc']] <- as.numeric(marker.tab[['log2fc']])*mult
            marker.tab[['test.stat']] <- as.numeric(marker.tab[['test.stat']])*mult

            marker.tab <- marker.tab[order(marker.tab$test.stat,decreasing=T),]
        }

        ## label.2 <- as.numeric(input$geneGroupSecond)+1

        ## label <- sprintf('%s.%s',label.1,label.2)

        ## marker.tab <- fromJSON(sprintf('http://localhost:3000/pairwise_markers/%s',label))

        ## rownames(marker.tab) <- marker.tab[,1]
        ## colnames(marker.tab) <- c('ens_id','p-val','unique-fcs','total-fcs')

        ## marker.tab <- cbind(cellmarker.info[rownames(marker.tab),],marker.tab[,c(2,3,4)],stringsAsFactors=F)

        ## marker.tab <- marker.tab[,c(-3)]
        ## colnames(marker.tab) <- c('Symbol','ID','gene-type','tissueType','cancerType','cellType','cellName','p-val','unique-fcs','total-fcs')

        marker.tab[['p-val']] <- as.numeric(marker.tab[['p-val']])

        return(marker.tab)

        ## if(input$geneGroup == 'all'){
        ##     dt <- marker.list[['all']]
        ## } else{
        ##     dt <- marker.list[[as.double(input$geneGroup)+1]]
        ## }

    })

    getSampleTable <- reactive({

        selectionRes <- selectionVec()
        louvainVec <- louvainVec()
        meshResults <- meshVec()
        tissueResults <- tissueVec()
        doidResults <- doidVec()
        efoResults <- efoVec()
        dt <- tsne.meta[,c(1,2,3,4,6,7,8,9)]
        colnames(dt) <- c('proj_id','samp_id','run_id','sample_type','url','tcga_case_id','tissue_general','tissue_detail')

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
        } else if(grepl('Louvain',input$sampleTableControl)){
            louvain.id <- gsub('Louvain Cluster ','',input$sampleTableControl)

            dt <- subset(dt,Louvain==as.numeric(louvain.id))
        }

        rownames(dt) <- c()

        ## print(head(dt))
        ## print(input$sampleTableControl)

        return(dt)

    })

    output$sampleTableDownload <- downloadHandler(
        filename = function() {
            paste('samples-',Sys.Date(),'.csv',sep='')
        },
        content = function(con) {
            write.csv(getSampleTable()[,c(-5)],con)
        }
    )

    output$geneTableDownload <- downloadHandler(
        filename = function() {
            paste('gene-table-',input$geneGroup,'-',Sys.Date(),'.csv',sep='')
        },
        content = function(con) {

            gene.tab <- getGeneTable()
            gene.tab <- gene.tab[,c(-8,-9,-10)]
            
            write.csv(gene.tab,con)
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

    ## output$metadataBar <- renderPlot({
    output$metadataBar <- renderPlotly({
        
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

        ## print(head(plot.dat))

        ## output.plot <- ggplot() +
        ##     geom_bar(data=plot.dat,aes(x=Label,y=number,fill=Label),stat='identity') +
        ##     theme(panel.background= element_blank(),
        ##           axis.text.x = element_text(size=14,angle=45,vjust=0.5),
        ##           axis.text.y = element_text(size=18),
        ##           legend.text = element_text(size=24),
        ##           axis.title.x = element_blank(),
        ##           axis.title.y = element_text(size=20),
        ##           legend.title = element_blank(),
        ##           legend.position='none',
        ##           plot.title=element_text(hjust=0.9,size=22)
        ##           ) +
        ##     scale_x_discrete(labels=x.labels,limits=names(x.labels)) +
        ##     ylab('') +
        ##     ggtitle(plot.title)

        output.plot <- plot_ly(data=plot.dat,
                               x=~Label,
                               y=~number,
                               color=~Label) %>%
            plotly::config(modeBarButtonsToRemove=c('hoverCompareCartesian','resetScale2d','hoverClosestCartesian','toggleSpikelines','zoomIn2d','zoomOut2d','autoScale2d','pan2d','zoom2d'))
        output.plot

    })
        

    violinData <- reactive({

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
            geneScale <- input$geneScale
        })

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
            xGroup <- xGroup[!is.na(xGroup)]
        } else{
            xGroup <- colVarResults
        }

        if(yFactors == 'projection'){

            yGroup <- -1 * projectionVec

            ylabel <- 'Projection'

            ## plot.title=paste('Projection of ',projection.list[sampleInputSelected,'name'],sep='')

        } else{

            yGroup <- geneVecResults
            gene.label <- names(gene.choices)[gene.num.vec[gene.id]]

            ## plot.title=paste('Expression of ',gene.label,sep='')
        }

        if(any(xGroup == 'All')){
            ## used.names <- names(yGroup)
            ## plot.dat <- data.frame(x=as.factor(xGroup),y=yGroup)
            plot.dat <- data.frame(x=xGroup,y=yGroup,stringsAsFactors=F)
            sizeVec <- rep(3,nrow(plot.dat))
            top.groups <- NULL
        } else{

            used.names <- intersect(names(xGroup),names(yGroup))

            plot.dat <- data.frame(x=as.character(xGroup[used.names]),y=yGroup[used.names],stringsAsFactors=F)

        } 
        
        return(plot.dat)

    })

    onevent('change','violinGroup',function(event){

        ## if(!(user.selections$violinGroupRendered)){
        ##     ## print('escaped violin change')
        ##     return()
        ## }

        if(user.selections$firstViolinRender){
            print('escape first render')
            user.selections$firstViolinRender <- FALSE
            user.selections$colorButtonActive <- FALSE
            return()
        }

        ## print('fire violin change')
        ## print(input$violinGroup)
        ## print(user.selections$currentViolinGroups)

        ## print(all(input$violinGroup == user.selections$currentViolinGroups))

        ## print(user.selections$colorButtonActive)

        col.butt <- input$colorButton

        if(col.butt == 0){
            return()
        }

        ## isolate({
        violinGroups <- input$violinGroup
        violinData <- violinData()
        ## })

        selectionRes <- selectionVec()

        if(user.selections$colorButtonActive){

            n.traces <- length(user.selections$currentViolinGroups)

            ## print(user.selections$currentViolinGroups)
            ## print(violinGroups)
            ## print(n.traces)

            color.ramp <- colorRampPalette(c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F"))(length(violinGroups))

            ## i <- 0
            
            set.seed(1)

            new.order <- violinGroups

            ## print(head(subset(selectionRes,group=='Selection: foo')))

            plotlyProxy('violin',session) %>%
                plotlyProxyInvoke(
                    'deleteTraces',
                    seq(0,n.traces-1)) %>%
                plotlyProxyInvoke(
                    'addTraces',
                    lapply(violinGroups,function(group.name) {
                        label <- gsub("(.{14,}?)\\s","\\1\n",group.name)
                        
                        ## i <<- i + 1

                        if(grepl('Selection',group.name)){

                            selectionSlice <- subset(selectionRes,group==group.name)

                            list(
                                frame=NULL,
                                x=rep(label,nrow(selectionRes)),
                                y=violinData[selectionSlice$samps,'y'],
                                xaxis='x',
                                yaxis='y',
                                type='violin',
                                showlegend=TRUE,
                                name=label)

                        } else if(group.name == 'Unlabelled'){

                            other.data <- subset(violinData,!(x %in% violinGroups))

                            list(
                                frame=NULL,
                                x=rep('Unlabelled',nrow(other.data)),
                                y=other.data[,'y'],
                                xaxis='x',
                                yaxis='y',
                                type='violin',
                                showlegend=TRUE,
                                name='Unlabelled')

                        } else{

                            dat.slice <- subset(violinData,x==group.name)

                            list(
                                ## fillcolor=color.ramp[i],
                                frame=NULL,
                                ## line=list(color=color.ramp[i]),
                                ## marker=list(color=color.ramp[i],
                                ##             line=list(
                                ##                 color=color.ramp[i])
                                ##             ),
                                x=rep(label,nrow(dat.slice)),
                                y=dat.slice[,'y'],
                                xaxis='x',
                                yaxis='y',
                                type='violin',
                                showlegend=TRUE,
                                name=label)

                        }
                    })) %>%
    ## plotlyProxyInvoke(
    ##     'addTraces',{
    ##         other.data <- subset(violinData,!(x %in% violinGroups))
    ##     list(
    ##         frame=NULL,
    ##         x=rep('Unlabelled',nrow(other.data)),
    ##         y=other.data[,'y'],
    ##         xaxis='x',
    ##         yaxis='y',
    ##         type='violin',
    ##         showlegend=TRUE,
    ##         name='Unlabelled')
    ##     }) %>%
    plotlyProxyInvoke(
        'relayout',
        list(showlegend=TRUE,
             'xaxis.autorange'=TRUE,
             
             ## 'xaxis.catqegoryarray'=violinGroups,
             ## 'xaxis.range'=c(0.5,length(violinGroups)+0.5),
             colorway=color.ramp)
        
    )
            
        } else{
            
            n.traces <- length(user.selections$currentViolinGroups)

            new.order <- user.selections$currentViolinGroups

            ### ADD ###
            if(n.traces < length(violinGroups)){
                
                new.group <- setdiff(violinGroups,user.selections$currentViolinGroups)
                
                print('add')
                print(new.group)

                if(new.group == 'Unlabelled'){

                    new.order <- c(new.order,'Unlabelled')

                    other.data <- subset(violinData,!(x %in% violinGroups))

                    plotlyProxy('violin',session) %>%
                        plotlyProxyInvoke(
                            'addTraces',{
                                list(
                                    frame=NULL,
                                    x=rep('Unlabelled',nrow(other.data)),
                                    y=other.data[,'y'],
                                    xaxis='x',
                                    yaxis='y',
                                    type='violin',
                                    showlegend=TRUE,
                                    name='Unlabelled')
                            })
                    
                } else if('Unlabelled' %in% user.selections$currentViolinGroups){

                    other.data <- subset(violinData,!(x %in% violinGroups))

                    if(grepl('Selection',new.group)){
                        data.slice <- violinData[subset(selectionRes,group==new.group)$samps,]
                        data.slice$x <- new.group

                    } else{
                        data.slice <- subset(violinData,x == new.group)
                    }

                    new.order <- new.order[-length(new.order)]
                    new.order <- c(new.order,new.group,'Unlabelled')
                    
                    plotlyProxy('violin',session) %>%
                        plotlyProxyInvoke('deleteTraces',c(n.traces-1)) %>%
                        plotlyProxyInvoke(
                            'addTraces',{
                                list(frame=NULL,
                                     x=rep(new.group,nrow(data.slice)),
                                     y=data.slice[,'y'],
                                     xaxis='x',
                                     yaxis='y',
                                     type='violin',
                                     showlegend=TRUE,
                                     name=new.group)
                            }) %>%
                            plotlyProxyInvoke(
                                'addTraces',{
                                list(frame=NULL,
                                     x=rep('Unlabelled',nrow(other.data)),
                                     y=other.data[,'y'],
                                    xaxis='x',
                                    yaxis='y',
                                    type='violin',
                                    showlegend=TRUE,
                                    name='Unlabelled')
                                })
                } else{

                    if(grepl('Selection',new.group)){
                        data.slice <- violinData[subset(selectionRes,group==new.group)$samps,]
                        data.slice$x <- new.group

                    } else{
                        data.slice <- subset(violinData,x == new.group)
                    }

                    new.order <- c(new.order,new.group)

                    plotlyProxy('violin',session) %>%                    
                        plotlyProxyInvoke(
                            'addTraces',{
                                list(frame=NULL,
                                     x=rep(new.group,nrow(data.slice)),
                                     y=data.slice[,'y'],
                                     xaxis='x',
                                     yaxis='y',
                                     type='violin',
                                     showlegend=TRUE,
                                     name=new.group)
                            })
                }

            } else{ ############# REMOVE

                new.group <- setdiff(user.selections$currentViolinGroups,violinGroups)

                print('remove')
                print(new.group)

                if(new.group == 'Unlabelled'){

                    new.order <- new.order[-n.traces]

                    plotlyProxy('violin',session) %>%
                        plotlyProxyInvoke('deleteTraces',n.traces-1)
                            
                } else if('Unlabelled' %in% user.selections$currentViolinGroups){

                    other.data <- subset(violinData,!(x %in% violinGroups))
                    remove.ind <- which(user.selections$currentViolinGroups == new.group)-1
                    ## data.slice <- subset(violinData,x == new.group)
                    
                    new.order <- new.order[-(remove.ind+1)]

                    plotlyProxy('violin',session) %>%
                        plotlyProxyInvoke('deleteTraces',c(remove.ind,(n.traces-1))) %>%
                        plotlyProxyInvoke(
                            'addTraces',{
                                list(frame=NULL,
                                     x=rep('Unlabelled',nrow(other.data)),
                                     y=other.data[,'y'],
                                     xaxis='x',
                                     yaxis='y',
                                     type='violin',
                                     showlegend=TRUE,
                                     name='Unlabelled')
                           })

                } else{

                    remove.ind <- which(user.selections$currentViolinGroups == new.group)-1
                    
                    new.order <- new.order[-(remove.ind+1)]

                    plotlyProxy('violin',session) %>%
                        plotlyProxyInvoke('deleteTraces',remove.ind) ## %>%
                        ## plotlyProxyInvoke(
                        ##     'relayout',
                        ##     list(showlegend=TRUE,
                        ##          'xaxis.categoryarray'=violinGroups,
                        ##          'xaxis.range'=c(-0.5,length(violinGroups)-0.5),
                        ##          colorway=color.ramp)
        
            ## )
                    
                }

            }

        }

        print(new.order)

        ## user.selections$currentViolinGroups <- violinGroups        
        user.selections$currentViolinGroups <- new.order
        
        user.selections$colorButtonActive <- FALSE

    })

    

    ## output$violin <- renderPlot({
    output$violin <- renderPlotly({
        ## input$colorButton

        isolate({

            col.butt <- input$colorButton
            violinData <- violinData()
            violinXGroup <- input$violinXFactors            

            selectionRes <- selectionVec()
            selectionRows <- input$selectionList_rows_selected
            selectionNames <- sapply(selectionRows,function(x) paste('Selection:',as.character(user.selections$selection.datalist[x,'name'])))

            user.selections$firstViolinRender <- TRUE

        })

        violinData <- subset(violinData, x != 'NA')

        if(any(violinData$x == 'All')){
            top.groups <- 'All'

            plot.dat <- violinData
            
        } else if(violinXGroup == 'Louvain'){
            group.means <- aggregate(violinData[,'y'],list(violinData[,'x']),function(x) mean(exp(x)))            
            group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]
            num.selections <- length(selectionNames)

            if(length(selectionNames) > 0){
                num.selections <- length(selectionNames)
            } else{
                ## selectionRes <- data.frame()
                num.selections <- 0
            }

            top.groups <- group.means[1:(10-num.selections)]

            top.groups <- c(top.groups,selectionNames,'Unlabelled')

            plot.dat <- violinData
            plot.dat[!(plot.dat$x %in% top.groups),'x'] <- 'Unlabelled'

        } else if(violinXGroup =='tissue_general'){

            top.groups <- unique(violinData$x)
            plot.dat <- violinData
            plot.dat[!(plot.dat$x %in% top.groups),'x'] <- 'Unlabelled'

        } else if(violinXGroup == 'tissue_detail'){

            if(length(selectionNames) > 0){
                num.selections <- length(selectionNames)
            } else{
                num.selections <- 0
            }

            group.means <- aggregate(violinData[,'y'],list(violinData[,'x']),function(x) mean(exp(x)))
            group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]
            group.means <- group.means[! group.means %in% c('Skin','Soft Tissue','Thymus')]

            top.groups <- group.means[1:(13 - num.selections)]

            plot.dat <- violinData
            plot.dat[!(plot.dat$x %in% top.groups),'x'] <- 'Unlabelled'

            top.groups <- c(top.groups,selectionNames,'Unlabelled')

        } else{

            labelled.dat <- subset(violinData,x != 'unlabelled')

            group.cnts <- table(labelled.dat[,'x'])
            used.groups <- names(which(group.cnts>(nrow(labelled.dat) / 100)))
            labelled.dat <- subset(labelled.dat,x %in% used.groups)

            if(length(selectionNames) > 0){
                num.selections <- length(selectionNames)
            } else{
                num.selections <- 0
            }

            if((length(unique(violinData$x)) + num.selections > 10)){

                group.means <- aggregate(labelled.dat[,'y'],list(labelled.dat[,'x']),mean)
                group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]

                top.groups <- group.means[1:(10 - num.selections)]

                plot.dat <- violinData
                plot.dat[!(plot.dat$x %in% top.groups),'x'] <- 'Unlabelled'

                top.groups <- c(top.groups,selectionNames,'Unlabelled')
            }
        }
        
        ##     used.names <- intersect(names(xGroup),names(yGroup))
        ##     plot.dat <- data.frame(x=as.character(xGroup[used.names]),y=yGroup[used.names],stringsAsFactors=F)

        ##     labelled.dat <- subset(plot.dat,x != 'unlabelled')

        ##     group.cnts <- table(labelled.dat[,'x'])
        ##     used.groups <- names(which(group.cnts>(nrow(labelled.dat) / 100)))
        ##     labelled.dat <- subset(labelled.dat,x %in% used.groups)


        ##     if(!is.null(selectionRes)){

        ##         selectionRes <- selectionRes[selectionRes$samps %in% used.names,]

        ##         num.selections <- length(unique(selectionRes$group))

        ##         ## plot.dat <- rbind(plot.dat,data.frame(x=as.factor(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))

        ##     } else{

        ##         selectionRes <- data.frame()
        ##         num.selections <- 0
        ##     }

        ##     sizeVec <- rep(3,nrow(plot.dat))

        ##     top.groups <- NULL

        ##     if((length(unique(xGroup)) + num.selections) > 10){

        ##         group.means <- aggregate(labelled.dat[,'y'],list(labelled.dat[,'x']),mean)
        ##         group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]

        ##         ## group.means <- group.means[group.means %in% used.groups]

        ##         top.groups <- group.means[1:(10 - num.selections)]

        ##         plot.dat[which(! (plot.dat$x %in% as.character(top.groups))),'x'] <- 'Other'

        ##         sizeVec[which(plot.dat$x != 'Other')] <- 5
                
        ##         top.groups <- c(top.groups,unique(selectionRes$group),'Other')

        ##     }

            

        ##     plot.dat <- rbind(plot.dat,data.frame(x=as.character(selectionRes$group),y=yGroup[as.character(selectionRes$samps)]))
            
        ##     sizeVec[which(plot.dat$x != 'Other' & plot.dat$x != 'unlabelled')] <- 5

        ## } 




        ## x.labels <- gsub("(.{14,}?)\\s","\\1\n",unique(plot.dat$x))
        ## names(x.labels) <- unique(plot.dat$x)

        ## if('Other' %in% x.labels){
        ##     x.labels <- x.labels[x.labels != 'Other']
        ##     x.labels <- c(x.labels,'Other'='Other')
        ## } else if('unlabelled' %in% x.labels){
        ##     x.labels <- x.labels[x.labels != 'unlabelled']
        ##     x.labels <- c(x.labels,'unlabelled'='Unlabelled')
        ## } else if('Eye' %in% x.labels){
        ##     x.labels <- x.labels[x.labels != 'Eye']
        ##     x.labels <- c(x.labels,'Eye'='Eye')
        ## }

        ## if(!is.null(top.groups)){

        ##     x.labels <- x.labels[top.groups]

        ## }
        
        ## min.y <- min(subset(plot.dat,y>0)$y)

        ## if(yTransform == 'log2gene'){
        ##     plot.dat$y <- log2(plot.dat$y)
        ##     ylabel <- sprintf('log2(%s)',ylabel)
        ## } else if(yTransform == 'log2gene1'){
        ##     plot.dat$y <- log2(plot.dat$y+1)
        ##     ylabel <- sprintf('log2(%s + 1)',ylabel)
        ## }

        ## if(yFactors == 'gene'){
        ##     if(geneScale == 'log2gene'){
        ##         ylabel <- 'log2(GENE)'
        ##         ## plot.dat$y <- log2(plot.dat$y)
        ##         plot.dat$y <- log2(2^(plot.dat$y)-1)
        ##     } else if(geneScale == 'log2gene1'){
        ##         ylabel <- 'log2(GENE + 1)'
        ##         ## plot.dat$y <- log2(plot.dat$y+1)
        ##         ## plot.dat$y <- plot.dat$y
        ##     } else{
        ##         ylabel <- 'GENE'
        ##         plot.day$y <- 2^(plot.dat$y) - 1
        ##     }
        ## }

        ## output.plot <- ggplot() +
        ##     ## geom_violin(data=plot.dat,aes(x=x,y=log2(abs(y)+0),fill=x),scale='width') +
        ##     ## geom_jitter(data=plot.dat,aes(x=x,y=log2(abs(y)+0)),width=0.1,alpha=0.1,colour='gray',size=sizeVec) +
        ##     geom_violin(data=plot.dat,aes(x=x,y=y+0,fill=x),scale='width')+
        ##     geom_jitter(data=plot.dat,aes(x=x,y=y+0),width=0.1,alpha=0.1,colour='gray',size=sizeVec) +
        ##     theme(##panel.background= element_blank(),
        ##         axis.text.x = element_text(size=14,angle=45,vjust=0.5),
        ##         axis.text.y = element_text(size=18),
        ##         axis.title.x = element_blank(),
        ##         axis.title.y = element_text(size=20),
        ##         legend.title = element_blank(),
        ##         legend.position='none',
        ##         plot.title=element_text(hjust=1.0,size=22)
        ##           ) +
        ##     ylab(ylabel) +
        ##     ggtitle(plot.title) +
        ##     scale_x_discrete(labels=x.labels,limits=names(x.labels))

            ## scale_
            ## geom_violin(data=dataResults, x = ~y1, y = ~y2,mode="markers",type='scatter',color = colVarResults,text=colVarResults)
            ## config(p = .,modeBarButtonsToRemove = c("zoom2d",'toImage','autoScale2d','hoverClosestGl2d'),collaborate=FALSE,cloud=FALSE) %>%
            ## config(collaborate=FALSE) %>%

        updateSelectizeInput(session,
                             'violinGroup',
                             choices=c(unique(violinData$x),unique(selectionRes$group),'Unlabelled'),
                             selected=top.groups,
                             label=NULL)

        isolate({
            user.selections$currentViolinGroups <- unlist(top.groups)
            ## user.selections$violinGroupRendered <- TRUE
        })

        ## print('violllin')

        print('plotly render violin')

        output.plot <- plot_ly() %>%
            plotly::config(modeBarButtonsToRemove=c('hoverCompareCartesian','resetScale2d','hoverClosestCartesian','toggleSpikelines','zoomIn2d','zoomOut2d','autoScale2d','pan2d','zoom2d'))

        ## print(unlist(top.groups))

        for(group in unlist(top.groups)){

            if(!grepl('Selection',group)){
                output.plot <- output.plot %>%
                    add_trace(data=subset(plot.dat,x==group),
                              x=~x,
                              y=~y,
                              type='violin',
                              name=group,
                              hoveron='violin')#,
                              ## color=~y)
            } else{

                select.dat <- violinData[subset(selectionRes,group==group)$samps,]
                select.dat$x <- group

                output.plot <- output.plot %>%
                    ## add_trace(data=violinData[subset(selectionRes,group==group)$samps,],
                    add_trace(data=select.dat,
                              x=~x,
                              y=~y,
                              type='violin',
                              name=group,
                              hoveron='violin')#,
                              ## color=~y)
            }
        }

        ## output.plot <- plot_ly() %>%
        ##     add_trace(

        ## output.plot <- plot_ly(plot.dat, x=~x,y=~y,
        ##                        color = ~x,
        ##                        type='violin',
        ##                        source='violin') %>%
        ##     plotly::config(modeBarButtonsToRemove=c('hoverCompareCartesian','resetScale2d','hoverClosestCartesian','toggleSpikelines','zoomIn2d','zoomOut2d','autoScale2d','pan2d','zoom2d'))

        output.plot <- output.plot %>%
            plotly::layout(
                'xaxis' = list(
                    'title'=list(text=''),
                    showticklabels=FALSE,
                    categoryorder='trace'))
                    ## categoryarray=unlist(top.groups)))

        output.plot

    })

    ## observe({

    ##     a <- input$colorButton

    ##     if(a > 0){

    ##         print('foo')
    ##         plotlyProxy('tsne',session) %>%
    ##             plotlyProxyInvoke(
    ##                 'restyle',
    ##                 list('marker.color'=list(sample(c('red','blue'),34759,replace=T)),
    ##                      'marker.showscale'=TRUE))

            
    ##                 ## list('marker.color'=list(runif(n=34759)),
    ##                 ##      'marker.colorscale'='Viridis',
    ##                 ##      'marker.showscale'=TRUE))

            
    ##     }
        
    ## })
    
        ## plotlyProxy('tsne',session) %>%
        ##     plotlyProxyInvoke(
        ##         'restyle',
        ##         list('marker.color'='red'))})
    
    observe({

        col.butt <- input$colorButton

        if(col.butt == 0){
            return()
        }

        isolate({
            euclid.file <- input$euclid_input
            spear.file <- input$spearman_input
            Euclid.pca.file <- input$euclid_pca_input
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
            geneScale <- input$geneScale
            layout.dat <- event_data('plotly_relayout',source='tsne')

            violinData <- violinData()
            violinXGroup <- input$violinXFactors
            user.selections$colorButtonActive <- TRUE
            selectionRes <- selectionVec()
            selectionRows <- input$selectionList_rows_selected
        })

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
            ## colVarPlot <- log2(geneVecResults+1)

            if(geneScale == 'log2gene'){
                ## ylabel <- 'log2(GENE)'
                geneVecResults <- log2(geneVecResults)
                geneVecResults <- log2(2^geneVecResults -1 )
            } else if(geneScale == 'log2gene1'){
                ## ylabel <- 'log2(GENE + 1)'
                ## geneVecResults <- log2(geneVecResults + 1)
                geneVecResults <- geneVecResults
            }
            else{
                geneVecResults <- 2^geneVecResults - 1
                ## ylabel <- 'GENE'
            }

            inf.gene.vals <- is.infinite(geneVecResults)
            inf.gene.sum <- sum(inf.gene.vals)

            if(inf.gene.sum > 0){
                showNotification(sprintf('Removing %d points with -Inf gene value',inf.gene.sum),type='error',duration=10)
            }
            geneVecResults <- geneVecResults[!inf.gene.vals]
            colVarPlot <- geneVecResults
            ## colVarPlot <- colVarPlot[!is.infinite(colVarPlot)]
            dataResults <- dataResults[names(colVarPlot),]

            

        } else if(colorFactors == 'KMeans'){
            colVarPlot <- as.factor(kmeansVec[rownames(dataResults)])

        } else if(colorFactors == 'Louvain'){
            ## colVarPlot <- as.factor(louvainVec[rownames(dataResults)])
            colVarPlot <- louvainVec[rownames(dataResults)]
            colVarPlot[is.na(colVarPlot)] <- -1
            colVarPlot <- as.factor(colVarPlot)

        } else if(colorFactors == 'Projection'){
            ## colVarPlot <- log2(projectionVec+1)
            colVarPlot <- projectionVec
            ## colVarPlot <- -1*projectionVec

            ## bottom.ramp <- colorRampPalette(c('#D73027','#FFFFBF'))
            ## bottom.ramp <- colorRampPalette(c('#FF0000','#FFFFBF'))
            bottom.ramp <- colorRampPalette(c('#152238','#3792CB'))
            top.ramp <- colorRampPalette(c('#3792CB','#FF0000'))

            color.ramp <- c(bottom.ramp(50),top.ramp(50))
            
            ## bottom.ramp <- colorRampPalette(c('#FF0000','#ADD8E6'))

            ## top.ramp <- colorRampPalette(c('#FFFFBF','#4575B4'))
            ## top.ramp <- colorRampPalette(c('#FFFFBF','#FFFFF8'))
            ## top.ramp <- colorRampPalette(c('#ADD8E6','#F6FBFC'))

            ## min.point <- min(min(colVarPlot),0.1439294) ### using true minimum, is waaaaaay too low
            ## min.point <- min(min(colVarPlot),13.51984) # 1st percentile
            ## ## min.point <- min(min(colVarPlot),10.85842) # testing, may be too small
            ## ## max.point <- max(max(colVarPlot),25.65136)
            ## max.point <- 25.65136
            
            ## mean.point <- min(mean(colVarPlot),16.69396)
            ## mean.point <- 16.69396
            ## mean.point <- mean(colVarPlot)

            ## print(min.point)
            ## print(max.point)

            ## min.point <- min(colVarPlot)
            ## max.point <- max(colVarPlot)

            ## mid.point <- round((mean.point - min.point) / (max.point - min.point)*100)

            ## print(mid.point)

            ## mid.point <- round((mean(colVarPlot) - min(colVarPlot)) / (max(colVarPlot) - min(colVarPlot)) * 100)

            ## color.ramp <- c(bottom.ramp(mid.point),top.ramp(100-mid.point))

            ## sample.min.point <- round((min(colVarPlot) - min.point) / (max.point - min.point)*100)

            ## print(sample.min.point)

            ## color.ramp <- color.ramp[sample.min.point:100]
            
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
            
            color.ramp <- colorRampPalette(c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F"))(length(unique(colVarPlot))) #,
                                   ## colorRampPalette(brewer.pal(12,"Set3"))(length(unique(colVarPlot))),
                                   ## sample(unique(colVarPlot)
            ## color.ramp['unlabelled'] <- '#F0F0F0'
            
        }

        if(!is.null(colVarPlot) & colorFactors != 'Gene' & colorFactors != 'Projection'){
            
            i <- 1

            ## isolate({user.selections$tsne.traces <- list()})

            tsne.traces <- list()

            for(label in sort(unique(colVarPlot))){

                inds <- which(colVarPlot == label)

                tsne.traces[[i]] <- rownames(dataResults)[inds]
                
                color.ramp[i] <- if(label=='unlabelled') '#F0F0F0' else color.ramp[i]
        
                ## print(head(user.selections$tsne.traces[[i]]

                i <- i+1
            }


            
            ## plotlyProxy('tsne',session) %>%
            ##     plotlyProxyInvoke(
            ##         'restyle',
            ##         'marker.color',runif(length(colVarPlot)))
            



            isolate({
                n.traces <- length(user.selections$tsne.traces)

                ## color.ramp <- colorRampPalette(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"))(n.traces)

                plotlyProxy('tsne',session) %>%
                    plotlyProxyInvoke(
                        'deleteTraces',
                        seq(0,n.traces-1)) %>%
                    plotlyProxyInvoke(
                        'addTraces',
                        lapply(tsne.traces,function(x) {
                            
                            name <- unname(colVarPlot[x[1]])
                            
                            list(x=dataResults[x,'y1'],
                                y=dataResults[x,'y2'],
                                mode='markers',
                                type='scattergl',
                                showlegend=(name!='unlabelled'),
                                opacity=if(name!='unlabelled') 0.8 else 0.25,
                                name=name,
                                marker=list(size=markerSize),
                                hoverinfo='text',
                                hovertext=paste(dataResults[x,'hover.labels'],
                                                '<br>Gene: ',dataResults[x,'Gene'],
                                                '<br>Annotation: ',colVarAnnotation[x]))
                        }))  %>%
                    plotlyProxyInvoke(
                        'relayout',
                        list(showlegend=TRUE,
                             colorway=color.ramp))
                ## colorway='Set3'))
                ##         colorway= c('#f3cec9', '#e7a4b6', '#cd7eaf', '#a262a9', '#6f4d96', '#3d3b72', '#182844')))
                
                
                            
                
                user.selections$tsne.traces <- tsne.traces
            })

        }
        else{
            isolate({

                n.traces <- length(user.selections$tsne.traces)

                tsne.traces <- list(rownames(dataResults))

                if(n.traces > 1){

                    plotlyProxy('tsne',session) %>%
                        plotlyProxyInvoke(
                            'deleteTraces',
                            seq(0,n.traces-1)) %>%
                        plotlyProxyInvoke(
                            'addTraces',
                            lapply(tsne.traces,function(x) {
                                
                                ## name <- unname(colVarPlot[x[1]])
                                
                                list(x=dataResults[x,'y1'],
                                     y=dataResults[x,'y2'],
                                     mode='markers',
                                     type='scattergl',
                                     showlegend=FALSE,
                                     ## opacity=if(name!='unlabelled') 0.8 else 0.25,
                                     ## name=name,
                                     marker=list(size=markerSize,
                                                 color=unname(colVarPlot),
                                                 colorscale=color.ramp,
                                                 cmin=min(colVarPlot),
                                                 cmax=max(colVarPlot),
                                                 showscale=TRUE),
                                     hoverinfo='text',
                                     hovertext=paste(dataResults[x,'hover.labels'],
                                                     '<br>Gene: ',dataResults[x,'Gene'],
                                                     '<br>Annotation: ',colVarAnnotation[x]))
                            }))  ## %>%
                        ## plotlyProxyInvoke(
                        ##     'relayout',
                        ##     list(showlegend=TRUE,
                        ##          colorway=NULL))
                }
                else{
                    
                    plotlyProxy('tsne',session) %>%
                        plotlyProxyInvoke(
                            'restyle',
                            list('marker.color'=list(unname(colVarPlot)),
                                 'marker.colorscale'=color.ramp,
                                 'marker.showscale'=TRUE))

                }

                user.selections$tsne.traces <- tsne.traces
            })
        }

############################
#### handle new violin plot
############################

        if(user.selections$violinGroupRendered){

            ## print('color button violin')

            selectionNames <- sapply(selectionRows,function(x) paste('Selection:',as.character(user.selections$selection.datalist[x,'name'])))

            violinData <- subset(violinData, x != 'NA')

            if(any(violinData$x == 'All')){
                top.groups <- 'All'

            } else if(violinXGroup == 'Louvain'){

                group.means <- aggregate(violinData[,'y'],list(violinData[,'x']),function(x) mean(exp(x)))

                group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]

                num.selections <- length(selectionNames)

                if(length(selectionNames) > 0){
                    num.selections <- length(selectionNames)
                } else{
                    ## selectionRes <- data.frame()
                    num.selections <- 0
                }

                top.groups <- group.means[1:(10-num.selections)]

                top.groups <- c(top.groups,selectionNames,'Unlabelled')

            } else if(violinXGroup=='tissue_general'){

                top.groups <- unique(violinData$x)

            } else if(violinXGroup=='tissue_detail'){

                if(length(selectionNames) > 0){
                    num.selections <- length(selectionNames)
                } else{
                    num.selections <- 0
                }

                group.means <- aggregate(violinData[,'y'],list(violinData[,'x']),function(x) mean(exp(x)))
                group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]
                group.means <- group.means[! group.means %in% c('Skin','Soft Tissue','Thymus')]

                top.groups <- group.means[1:(13 - num.selections)]

                ## print(selectionNames)

                top.groups <- c(top.groups,selectionNames,'Unlabelled')
            } else{

                labelled.dat <- subset(violinData,x != 'unlabelled')

                group.cnts <- table(labelled.dat[,'x'])
                used.groups <- names(which(group.cnts>(nrow(labelled.dat) / 100)))
                labelled.dat <- subset(labelled.dat,x %in% used.groups)

                if(length(selectionNames) > 0){
                    num.selections <- length(selectionNames)
                } else{
                    num.selections <- 0
                }

                if((length(unique(violinData$x)) + num.selections > 10)){

                    group.means <- aggregate(labelled.dat[,'y'],list(labelled.dat[,'x']),mean)
                    group.means <- group.means[['Group.1']][order(group.means$x,decreasing=TRUE)]

                    top.groups <- group.means[1:(10 - num.selections)]

                    top.groups <- c(top.groups,selectionNames,'Unlabelled')
                } else{

                    top.groups <- c(unique(labelled.dat[,'x']),selectionNames,'Unlabelled')

                }
            }

            ## print('top groups')
            ## print(top.groups)
            ## print(c(unique(violinData$x),unique(selectionRes$group),'Unlabelled'))
            

            updateSelectizeInput(session,
                                 'violinGroup',
                                 choices = c(unique(violinData$x),unique(selectionRes$group),'Unlabelled'),
                                 selected = top.groups,
                                 label=NULL)
        }##  else{
        ##     print('escaped color button')
        ## }
    })        

    output$tsne <- renderPlotly({

        ## input$colorButton

        isolate({
            euclid.file <- input$euclid_input
            spear.file <- input$spearman_input
            Euclid.pca.file <- input$euclid_pca_input
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
            geneScale <- input$geneScale
            layout.dat <- event_data('plotly_relayout',source='tsne')
            colorButton <- input$colorButton
        })


        
        ## print(colorFactors)

        ## color.ramp <- c('#800080','#FFFF00')

        color.ramp <- NULL

        ## if(input$colorButton == 0 | colorFactors == 'No Coloring'){
        if(colorButton == 0 | colorFactors == 'No Coloring'){
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
            ## colVarPlot <- log2(geneVecResults+1)

            if(geneScale == 'log2gene'){
                ## ylabel <- 'log2(GENE)'
                geneVecResults <- log2(geneVecResults)
                geneVecResults <- log2(2^geneVecResults -1 )
            } else if(geneScale == 'log2gene1'){
                ## ylabel <- 'log2(GENE + 1)'
                ## geneVecResults <- log2(geneVecResults + 1)
                geneVecResults <- geneVecResults
            }
            else{
                geneVecResults <- 2^geneVecResults - 1
                ## ylabel <- 'GENE'
            }

            inf.gene.vals <- is.infinite(geneVecResults)
            inf.gene.sum <- sum(inf.gene.vals)

            if(inf.gene.sum > 0){
                showNotification(sprintf('Removing %d points with -Inf gene value',inf.gene.sum),type='error',duration=10)
            }
            geneVecResults <- geneVecResults[!inf.gene.vals]
            colVarPlot <- geneVecResults
            ## colVarPlot <- colVarPlot[!is.infinite(colVarPlot)]
            dataResults <- dataResults[names(colVarPlot),]

            

        } else if(colorFactors == 'KMeans'){
            colVarPlot <- as.factor(kmeansVec[rownames(dataResults)])

        } else if(colorFactors == 'Louvain'){
            ## colVarPlot <- as.factor(louvainVec[rownames(dataResults)])
            colVarPlot <- louvainVec[rownames(dataResults)]
            colVarPlot[is.na(colVarPlot)] <- -1
            colVarPlot <- as.factor(colVarPlot)

        } else if(colorFactors == 'Projection'){
            ## colVarPlot <- log2(projectionVec+1)
            colVarPlot <- projectionVec
            ## colVarPlot <- -1*projectionVec

            ## bottom.ramp <- colorRampPalette(c('#D73027','#FFFFBF'))
            ## bottom.ramp <- colorRampPalette(c('#FF0000','#FFFFBF'))
            bottom.ramp <- colorRampPalette(c('#152238','#3792CB'))
            top.ramp <- colorRampPalette(c('#3792CB','#FF0000'))

            color.ramp <- c(bottom.ramp(50),top.ramp(50))
            
            ## bottom.ramp <- colorRampPalette(c('#FF0000','#ADD8E6'))

            ## top.ramp <- colorRampPalette(c('#FFFFBF','#4575B4'))
            ## top.ramp <- colorRampPalette(c('#FFFFBF','#FFFFF8'))
            ## top.ramp <- colorRampPalette(c('#ADD8E6','#F6FBFC'))

            ## min.point <- min(min(colVarPlot),0.1439294) ### using true minimum, is waaaaaay too low
            ## min.point <- min(min(colVarPlot),13.51984) # 1st percentile
            ## ## min.point <- min(min(colVarPlot),10.85842) # testing, may be too small
            ## ## max.point <- max(max(colVarPlot),25.65136)
            ## max.point <- 25.65136
            
            ## mean.point <- min(mean(colVarPlot),16.69396)
            ## mean.point <- 16.69396
            ## mean.point <- mean(colVarPlot)

            ## print(min.point)
            ## print(max.point)

            ## min.point <- min(colVarPlot)
            ## max.point <- max(colVarPlot)

            ## mid.point <- round((mean.point - min.point) / (max.point - min.point)*100)

            ## print(mid.point)

            ## mid.point <- round((mean(colVarPlot) - min(colVarPlot)) / (max(colVarPlot) - min(colVarPlot)) * 100)

            ## color.ramp <- c(bottom.ramp(mid.point),top.ramp(100-mid.point))

            ## sample.min.point <- round((min(colVarPlot) - min.point) / (max.point - min.point)*100)

            ## print(sample.min.point)

            ## color.ramp <- color.ramp[sample.min.point:100]
            
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

        if(is.null(layout.dat)){
            
            ax.x <- ax.y <- list(
                title = "",
                zeroline = FALSE,
                showline = FALSE,
                showticklabels = FALSE,
                showgrid = FALSE,
                range = c(-110,110)
            )
            
        }
        else{

            ax.x <- list(
                title = "",
                zeroline = FALSE,
                showline = FALSE,
                showticklabels = FALSE,
                showgrid = FALSE,
                range=c(layout.dat[['xaxis.range[0]']],layout.dat[['xaxis.range[1]']]))

            ax.y <- list(
                title = "",
                zeroline = FALSE,
                showline = FALSE,
                showticklabels = FALSE,
                showgrid = FALSE,
                range=c(layout.dat[['yaxis.range[0]']],layout.dat[['yaxis.range[1]']]))

        }


        ## if(colorButton == 0 ){

        output.plot <- plot_ly(dataResults, x = ~y1, y = ~y2,mode="markers",type='scattergl',color = colVarPlot,
                               colors=color.ramp,
                               hoverinfo='text',
                               ## opacity=plot.opacity,
                               ## group = colVarPlot,
                               ## text=~paste('ID: ',run_id,
                               ##             '<br>Project: ',proj_id,
                               ##             '<br>Tissue General: ',tissue_general,
                               ##             '<br>Tissue Detail: ',tissue_detail,
                               ##             '<br>Sample Type: ',sample_type,
                               ##             '<br>Sample ID: ',samp_id,
                               ##             ## '<br>Kmeans: ',KMeans,
                               ##             '<br>Louvain: ',Louvain,
                               ##             '<br>Gene: ',Gene,
                               ##             '<br>Annotation: ',colVarAnnotation),
                               text=~paste(hover.labels,
                                           '<br>Gene: ',Gene,
                                           '<br>Annotation: ',colVarAnnotation),
                               marker = list(size=markerSize),
                               source='tsne') %>%
            layout(dragmode = "pan",
                   xaxis=ax.x,
                   yaxis=ax.y,
                   legend=list(font=list(family='sans-serif',size=11))) %>%
            ## modebar=list(modeBarButtonsToRemove=c('zoomIn2d','lasso2d')))  %>%
            ## plotly::config(modeBarButtonsToRemove=c('zoom2d'))
            plotly::config(modeBarButtonsToRemove=c('hoverCompareCartesian','resetScale2d','hoverClosestCartesian','toggleSpikelines'))

        ## config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d"))            

            ## config(output.plot,modeBarButtonsToRemove = c('zoom2d'))


            ## output.plot <- config(output.plot,modeBarButtonsToRemove = c('zoom2d'))

            ## config(p = .,modeBarButtonsToRemove = c("zoom2d",'toImage','autoScale2d','hoverClosestGl2d'),collaborate=FALSE,cloud=FALSE) %>%
            ## config(collaborate=FALSE) %>%


        output.plot
        ## }

    })
  
})
