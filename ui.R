library(shiny)
library(plotly)
library(shinyTree)
library(DT)
library(shinyjs)

meta.choices <- c('No Coloring' = 'No Coloring',
                  'Project ID' = 'proj_id',
                  'GTEx/TCGA Gross Tissue' = 'tissue_general',
                  'GTEx/TCGA Detailed Tissue' = 'tissue_detail',
                  'Sample Type' = 'sample_type',
                  'Mesh: Anatomy' = 'Mesh',
                  'Tissue' = 'Tissue',
                  'DOID' = 'DOID',
                  'EFO: Cultured Cells'='efo',
                  'Gene' = 'Gene',
                  ## 'Kmeans' = 'KMeans',
                  'Louvain' = 'Louvain',
                  'Projection' = 'Projection')




## louvain.vec <- readRDS('data/leiden_r25e-3_over50_pc3sd_mrmnorm_k40_sim_nosingles.RDS')
## louvain.vec <- readRDS('data/leiden_r25e-3_over50_pc3sd_poscounts_k40_sim_nosingles.RDS')
## louvain.vec <- readRDS('data/leiden_pca_r5e-3_pc3sd_tpm_log_k20_sim_nosingles.RDS')
louvain.vec <- readRDS('data/leiden_r9e-3_over50_pc3sd_tpm_log_k30_sim_90th_var_genes.RDS')
naked.louvain.choices <- sort(unique(louvain.vec))
louvain.choices <- sapply(naked.louvain.choices,function(x) sprintf('Louvain Cluster %s',x))
names(louvain.choices) <- louvain.choices

## wgcna.choices <- readRDS('data/wgcna_app_choices_pos_log_pc3sd_min100_ds2.RDS')

fluidPage(
    includeCSS('styles.css'),
    includeScript('jssrc/panel_click.js'),
    navbarPage('Mercator',id='topNavBar',
               tabPanel('t-SNE',value='tsnePanel',
                        absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=75,left=0,right='auto',bottom='auto',width='100%',height='95%',value='tsnePanel',
                                      plotlyOutput('tsne',width='100%',height='100%')
                                      )
                        ),
               tabPanel('Violin',value='violinPanel',
                                      absolutePanel(id='results',fixed=FALSE,draggable=FALSE,top=70,left=25,right='auto',bottom='auto',width='97%',height='75%',value='violinPanel',
                                                    plotlyOutput('violin',width='100%',height='100%'),

                                                    fluidRow(
                                                        column(11,offset=0,
                                                               selectizeInput('violinGroup',
                                                                              label=NULL,
                                                                              choices=NULL,
                                                                              multiple=TRUE,
                                                                              width='100%',
                                                                              options=list(
                                                                                  dropdownDirection='up')
                                                                              )
                                                               ))
                                      )
                        
                        ## absolutePanel(id='violinSelect',fixed=TRUE,draggable=FALSE,top=400,left=25,width='90%',bottom='auto',
                        ##               selectizeInput('violinGroup',
                        ##                              label=NULL,
                        ##                              choices=NULL,
                        ##                              multiple=TRUE,
                        ##                              width='100%'))

                        ),
               tabPanel('Bar Plot',value='barPanel',
                        absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=75,left=15,right='auto',bottom='auto',width='99%',height='85%',value='barPanel',
                                      plotlyOutput('metadataBar',width='100%',height='100%')
                                      )
                        ),

               tabPanel('Cluster Markers',value='genePanel',
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=72,left=435,right='auto',bottom='auto',width='auto',height='auto',
                                      tags$label(class='control-label', 'Cluster #s for marker table')
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=42,left=640,right='auto',bottom='auto',width='auto',height='auto',
                                      selectInput('geneGroup',
                                                  label='',
                                                  ## label='Cluster # for marker table',
                                                  ## choices = c('All'='all',naked.louvain.choices),
                                                  choices = naked.louvain.choices,
                                                  multiple=FALSE,
                                                  selected='All',
                                                  width='100px')
                                      ),

                        absolutePanel(fixed=FALSE,draggable=FALSE,top=72,left=750,right='auto',bottom='auto',width='auto',height='auto',
                                      tags$label(class='control-label', 'vs')
                                      ),

                        absolutePanel(fixed=FALSE,draggable=FALSE,top=42,left=780,right='auto',bottom='auto',width='auto',height='auto',
                                      selectInput('geneGroupSecond',
                                                  label='',
                                                  ## label='Cluster # for marker table',
                                                  choices = c('All'='all',naked.louvain.choices),
                                                  multiple=FALSE,
                                                  selected='All',
                                                  width='100px')
                                      ),
                                      
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=62,left=930,right='auto',bottom='auto',width='auto',height='auto',
                                      downloadButton('geneTableDownload',
                                                     label = 'Download gene table'
                                                     )
                                      ),
                        absolutePanel(id='results',fixed=FALSE,draggable=FALSE,top=120,left=50,right='auto',bottom='auto',width='99%',height='auto',value='genePanel',

                                    ## DT::DTOutput('geneTable',width='300px')
                                    DT::DTOutput('geneTable',width='90%')
                                    )
                        ),

               tabPanel('Sample Table',value='samplePanel',
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=72,left=435,right='auto',bottom='auto',width='auto',height='auto',
                                      tags$label(class='control-label', 'Subset data')
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=42,left=550,right='auto',bottom='auto',width='auto',height='auto',
                                      selectInput('sampleTableControl',
                                                  label='',
                                                  ## label='Cluster # for marker table',
                                                  choices = c('All'='All','Selection'='Selection',louvain.choices),
                                                  multiple=FALSE,
                                                  selected=0)
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=62,left=875,right='auto',bottom='auto',width='auto',height='auto',
                                      downloadButton('sampleTableDownload',
                                                     label = 'Download sample table'
                                                     )
                                      ),
                        absolutePanel(id='results',fixed=FALSE,draggable=FALSE,top=120,left=50,right='auto',bottom='auto',width='99%',height='auto',
                                      DT::DTOutput('sampleTable',width='90%')
                                      )
                        
                        )## ,
               ## tabPanel('Gene Clustering',value='wgcnaPanel',

               ##          absolutePanel(fixed=FALSE,draggable=FALSE,top=72,left=435, right='auto', bottom='auto',width='auto',height='auto',
               ##                        tags$label(class='control-label','Gene Cluster #')
               ##                        ),
               ##          absolutePanel(fixed=FALSE,draggable=FALSE,top=42,left=550,right='auto',bottom='auto',width='auto',height='auto',
               ##                        selectInput('wgcnaTabControl',
               ##                                    label='',
               ##                                    choices = c(wgcna.choices),
               ##                                    multiple=FALSE,
               ##                                    selected=1)
               ##                        ),
               ##          absolutePanel(id='results',fixed=FALSE,draggable=FALSE,top=120,left=50,right='auto',bottom='auto',width=600,height='auto',
               ##                        DT::DTOutput('wgcnaGeneTable',width='90%')
               ##                        ),
               ##          absolutePanel(id='results',fixed=FALSE,draggable=FALSE,top=120,left=675,right='auto',bottom='auto',width=600,height='auto',
               ##                        DT::DTOutput('wgcnaGoTable',width='90%')
               ##                        )

                        ## )
               ),
    
    absolutePanel(id='controlPanel', fixed=FALSE, draggable=TRUE, top=60, left=20, right='auto', bottom='auto', width=400, height='auto',style="z-index:100;",
                  
                  ## wellPanel(

                  tags$div(class='panel panel-primary',
                           tags$div(class='panel-heading',id='controlPanelHeading',
                                    
                                    useShinyjs(),
                                    h3("Plot Controls",class='panel-title'),
                                    tags$span(class='pull-right clickable',id='controlPanelClickable',tags$i(class='glyphicon glyphicon-chevron-up'))
                                    ),
                           tags$div(class='panel-body',

                                    ## actionButton('colorButton','Redraw Plot'),
                                    actionButton('colorButton',tags$h4('Redraw Plots')),


                                    tags$hr(class='section-divider'),

                                    ## tags$br(),
                                    ## tags$br(),

                                    uiOutput('plotlyClick'),
                                    
                                    ## tags$br(),
                                    tags$hr(class='section-divider'),

                                    selectizeInput('whichGene',
                                                   label='Gene selection',
                                                   choices=NULL),

                                    selectInput('geneScale',
                                                label='Gene Value Transformation',
                                                choices = c('log2(GENE)'='log2gene','log2(GENE+1)'='log2gene1','None'='none'),
                                                multiple=FALSE,
                                                selected=c('log2(GENE+1)'='log2gene1')
                                                ),

                                    tags$hr(class='section-divider'),



                                    ## tags$br(),

                                    selectInput('colorfactors',
                                                label = 't-SNE Color',
                                                choices = meta.choices,
                                                multiple=FALSE
                                                ),

                                    sliderInput('markerSize','t-SNE Marker Size',min=1,max=12,value=4),

                                    tags$hr(class='section-divider'),

                                    selectInput('violinXFactors',
                                                label= 'Violin plot X-axis',
                                                choices = meta.choices,
                                                multiple=FALSE
                                                ),

                                    selectInput('violinYFactors',
                                                label='Violin plot Y-axis',
                                                choices = c('Projection'='projection','Gene'='gene'),
                                                multiple=FALSE,
                                                selected=c('Gene'='gene')
                                                ),


                                    tags$hr(class='section-divider'),

                                    selectInput('barPlotFactor',
                                                label='Bar plot subset',
                                                choices = c('All','Selections',louvain.choices),
                                                multiple=FALSE
                                                ),

                                    selectInput('barPlotXaxis',
                                                ## label='X axis for Bar Plot',
                                                label = 'Bar plot X-axis',
                                                choices=c('No Coloring' = 'No Coloring',
                                                          'Project ID' = 'proj_id',
                                                          'GTEx/TCGA Gross Tissue' = 'tissue_general',
                                                          'GTEx/TCGA Detailed Tissue' = 'tissue_detail',
                                                          'Sample Type' = 'sample_type',
                                                          'Mesh: Anatomy' = 'Mesh',
                                                          'Tissue' = 'Tissue',
                                                          'DOID' = 'DOID',
                                                          'EFO: Cultured Cells'='efo',
                                                          'Louvain' = 'Louvain'),
                                                multiple=FALSE
                                                ),

                                    tags$hr(class='section-divider'),
                                    
                                    tags$b('Ontology selection'),
                                    tags$br(),

                                    tags$p('Mesh: Anatomy'),
                                    shinyTree('meshTree',theme='proton'),
                                    tags$br(),

                                    tags$p('Tissue Ontology'),
                                    shinyTree('tissueTree',theme='proton'),
                                    tags$br(),

                                    tags$p('DOID'),
                                    shinyTree('doidTree',theme='proton'),
                                    tags$br(),

                                    tags$p('EFO: Cultured cell'),
                                    shinyTree('efoTree',theme='proton'),
                                    tags$br(),
                                    

                                    tags$hr(class='section-divider'),
                                    ## textInput('markerSearch','Search Marker Table'),
                                    
                                    textInput('selectionName','Name lasso/box selection',value=''),

                                    actionButton('saveSelection','Save Selection'),

                                    tags$br(),

                                    tags$br(),

                                    tags$b('Selection Table'),
                                    DT::DTOutput('selectionList',width='17%'),

                                    tags$hr(class='section-divider'),

                                    fileInput('gene.vec','Color by correlation to sample',
                                              accept = c(
                                                  'text/tsv',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.tsv')
                                              ),
                                    
                                    tags$b('Sample Input Table'),
                                    DT::DTOutput('sampleInputTable',width='300px')



                                    
                                    )
                           )
                  )
)
