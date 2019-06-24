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




louvain.vec <- readRDS('data/leiden_pc3sd_r25e-3_vec.RDS')
naked.louvain.choices <- sort(unique(louvain.vec))
louvain.choices <- sapply(naked.louvain.choices,function(x) sprintf('Louvain Cluster %s',x))
names(louvain.choices) <- louvain.choices

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
                        absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=70,left=25,right='auto',bottom='auto',width='97%',height='90%',value='violinPanel',
                                      plotOutput('violin',width='100%',height='100%')
                                      )
                        ),
               tabPanel('Bar Plot',value='barPanel',
                        absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=75,left=15,right='auto',bottom='auto',width='99%',height='85%',value='barPanel',
                                      plotOutput('metadataBar',width='100%',height='100%')
                                      )
                        ),

               tabPanel('Gene Table',value='genePanel',
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=72,left=435,right='auto',bottom='auto',width='auto',height='auto',
                                      tags$label(class='control-label', 'Cluster # for marker table')
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=42,left=620,right='auto',bottom='auto',width='auto',height='auto',
                                      selectInput('geneGroup',
                                                  label='',
                                                  ## label='Cluster # for marker table',
                                                  choices = c('All'='all',naked.louvain.choices),
                                                  multiple=FALSE,
                                                  selected='All')
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=62,left=975,right='auto',bottom='auto',width='auto',height='auto',
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
                        
                        )
               ## ),
               ),
    
    absolutePanel(id='controlPanel', fixed=FALSE, draggable=TRUE, top=60, left=20, right='auto', bottom='auto', width=400, height='auto',
                  
                  ## wellPanel(

                  tags$div(class='panel panel-primary',
                           tags$div(class='panel-heading',id='controlPanelHeading',
                                    
                                    useShinyjs(),                           
                                    h3("Plot Controls",class='panel-title'),
                                    tags$span(class='pull-right clickable',id='controlPanelClickable',tags$i(class='glyphicon glyphicon-chevron-up'))
                                    ),
                           tags$div(class='panel-body',

                                    actionButton('colorButton','Redraw Plot'),
                                    
                                    tags$br(),
                                    tags$br(),

                                    sliderInput('markerSize','Adjust Marker Size',min=1,max=12,value=4),


                                    uiOutput('plotlyClick'),
                                    
                                    tags$br(),

                                    selectizeInput('whichGene',
                                                   label='Select a gene',
                                                   choices=NULL),


                                    fileInput('gene.vec','Color by euclidian distance to sample',
                                              accept = c(
                                                  'text/tsv',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.tsv')
                                              ),

                                    DT::DTOutput('sampleInputTable',width='300px'),

                                    tags$br(),

                                    selectInput('colorfactors',
                                                label = 'What to color by?',
                                                choices = meta.choices,
                                                multiple=FALSE
                                                ),

                                    selectInput('violinXFactors',
                                                label= 'Violin X-axis',
                                                choices = meta.choices,
                                                multiple=FALSE
                                                ),

                                    selectInput('violinYFactors',
                                                label='Violin Y-axis',
                                                choices = c('Projection'='projection','Gene'='gene'),
                                                multiple=FALSE,
                                                selected=c('Gene'='gene')
                                                ),

                                    selectInput('barPlotFactor',
                                                label='Bar plot group',
                                                choices = c('All','Selections',louvain.choices),
                                                multiple=FALSE
                                                ),

                                    selectInput('barPlotXaxis',
                                                label='X axis for Bar Plot',
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

                                    shinyTree('meshTree',theme='proton'),

                                    shinyTree('tissueTree',theme='proton'),

                                    shinyTree('doidTree',theme='proton'),

                                    shinyTree('efoTree',theme='proton'),
                                    

                                    ## textInput('markerSearch','Search Marker Table'),
                                    
                                    textInput('selectionName','Selection Name',value=''),

                                    actionButton('saveSelection','Save Selection'),

                                    DT::DTOutput('selectionList',width='17%')
                                    
                                    )
                           )
                  )
)

